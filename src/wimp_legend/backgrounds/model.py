"""Generic observed-energy background counts; empirical inputs are external."""
from dataclasses import dataclass
from typing import Callable, Protocol
import numpy as np
from .._analysis import IntegrationPolicy, array, counts, edges, elementwise, output, sample, scalar, text


class BackgroundModel(Protocol):
    def expected_counts(self, observed_bin_edges_keV, exposure_kg_day):
        ...


@dataclass(frozen=True)
class ZeroBackground:
    def expected_counts(self, observed_bin_edges_keV, exposure_kg_day):
        bins = edges(observed_bin_edges_keV)
        scalar(exposure_kg_day, 'exposure_kg_day')
        return np.zeros(len(bins)-1)


@dataclass(frozen=True)
class FlatBackground:
    rate_per_kg_day_keV_obs: float

    def __post_init__(self):
        object.__setattr__(self, 'rate_per_kg_day_keV_obs', scalar(
            self.rate_per_kg_day_keV_obs, 'rate_per_kg_day_keV_obs'))

    def expected_counts(self, observed_bin_edges_keV, exposure_kg_day):
        bins = edges(observed_bin_edges_keV)
        exposure = scalar(exposure_kg_day, 'exposure_kg_day')
        with np.errstate(over='raise', invalid='raise'):
            return output(exposure * self.rate_per_kg_day_keV_obs * np.diff(bins))


@dataclass(frozen=True)
class CallableBackground:
    callback: Callable
    policy: IntegrationPolicy = IntegrationPolicy()
    breakpoints_keV_obs: tuple[float, ...] = ()
    provenance: str = 'user-supplied; unverified'

    def __post_init__(self):
        if not callable(self.callback) or not isinstance(self.policy, IntegrationPolicy):
            raise TypeError('background callback/policy invalid')
        points = array(tuple(self.breakpoints_keV_obs), 'breakpoints_keV_obs', nonnegative=False)
        if points.ndim != 1:
            raise ValueError('breakpoints_keV_obs must be a vector')
        object.__setattr__(self, 'breakpoints_keV_obs', tuple(points))
        object.__setattr__(self, 'provenance', text(self.provenance, 'provenance'))

    def differential_rate(self, E_obs_keV):
        return elementwise(self.callback, E_obs_keV, 'background rate', nonnegative_input=False)

    def expected_counts(self, observed_bin_edges_keV, exposure_kg_day):
        bins = edges(observed_bin_edges_keV)
        exposure = scalar(exposure_kg_day, 'exposure_kg_day')
        if exposure == 0:
            return np.zeros(len(bins)-1)
        return output([exposure*self.policy.integrate(
            lambda e: sample(self.callback, e, 'background rate'), lo, hi,
            breakpoints=self.breakpoints_keV_obs) for lo, hi in zip(bins[:-1], bins[1:])])


@dataclass(frozen=True)
class BinnedBackground:
    """Exact-edge binned input, with explicit exposure semantics.

    normalization='per_unit_exposure': values in events/(kg day), no reference
    exposure. normalization='fixed_expected_counts': values are events for the
    required exposure_kg_day, and a different requested exposure is rejected.
    """
    observed_bin_edges_keV: tuple[float, ...]
    values: tuple[float, ...]
    normalization: str
    exposure_kg_day: float | None = None
    provenance: str = 'user-supplied; unverified'

    def __post_init__(self):
        bins = edges(self.observed_bin_edges_keV)
        values = counts(self.values, 'background values')
        if values.shape != (len(bins)-1,):
            raise ValueError('background values must match observed bins')
        mode = text(self.normalization, 'normalization')
        if mode == 'per_unit_exposure':
            if self.exposure_kg_day is not None:
                raise ValueError('per-unit-exposure values must not specify a fixed exposure')
        elif mode == 'fixed_expected_counts':
            exposure = scalar(self.exposure_kg_day, 'exposure_kg_day')
            if exposure == 0 and np.any(values > 0):
                raise ValueError('nonzero fixed counts cannot represent zero exposure')
            object.__setattr__(self, 'exposure_kg_day', exposure)
        else:
            raise ValueError('unknown background normalization')
        object.__setattr__(self, 'observed_bin_edges_keV', tuple(bins))
        object.__setattr__(self, 'values', tuple(values))
        object.__setattr__(self, 'normalization', mode)
        object.__setattr__(self, 'provenance', text(self.provenance, 'provenance'))

    def expected_counts(self, observed_bin_edges_keV, exposure_kg_day):
        bins = edges(observed_bin_edges_keV)
        exposure = scalar(exposure_kg_day, 'exposure_kg_day')
        if not np.array_equal(bins, self.observed_bin_edges_keV):
            raise ValueError('binned background requires exact edges; no implicit rebinning')
        if self.normalization == 'fixed_expected_counts':
            if exposure != self.exposure_kg_day:
                raise ValueError('fixed counts belong to a different exposure')
            return np.array(self.values)
        return output(np.asarray(self.values)*exposure)


@dataclass(frozen=True)
class CompositeBackground:
    components: tuple[BackgroundModel, ...]

    def __post_init__(self):
        values = tuple(self.components)
        if any(not callable(getattr(component, 'expected_counts', None)) for component in values):
            raise TypeError('background components must implement expected_counts')
        object.__setattr__(self, 'components', values)

    def expected_counts(self, observed_bin_edges_keV, exposure_kg_day):
        bins = edges(observed_bin_edges_keV)
        exposure = scalar(exposure_kg_day, 'exposure_kg_day')
        result = np.zeros(len(bins)-1)
        for component in self.components:
            values = counts(component.expected_counts(bins.copy(), exposure), 'background component counts')
            if values.shape != result.shape:
                raise ValueError('background component has incompatible count shape')
            with np.errstate(over='raise', invalid='raise'):
                result += values
        return output(result)
