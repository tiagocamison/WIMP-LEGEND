"""Explicit truth/observed background inputs, coupled only to response ports."""
from dataclasses import dataclass
from typing import Callable, Protocol

import numpy as np

from .._analysis import array, counts, edges, output, scalar, text
from ..energy import EnergyDomain


class TruthResponse(Protocol):
    @property
    def input_energy_domain(self) -> EnergyDomain:
        ...

    def expected_counts(self, truth_spectrum, observed_bin_edges_keV, exposure_kg_day,
                        truth_energy_range_keV, *, truth_breakpoints=()):
        ...


@dataclass(frozen=True)
class MonoenergeticLine:
    """Integrated intensity in events/(kg day), with an explicit energy domain.

    Truth intensities precede detector efficiencies. Observed intensities already
    include all detector/selection effects. This is a point measure, not a PDF.
    """
    energy_keV: float
    rate_per_kg_day: float
    energy_domain: EnergyDomain

    def __post_init__(self):
        if not isinstance(self.energy_domain, EnergyDomain):
            raise TypeError('energy_domain must be an EnergyDomain')
        value = array(self.energy_keV, 'line energy_keV',
                      nonnegative=self.energy_domain is not EnergyDomain.OBSERVED)
        if value.ndim != 0:
            raise ValueError('line energy_keV must be a scalar')
        object.__setattr__(self, 'energy_keV', float(value))
        object.__setattr__(self, 'rate_per_kg_day', scalar(self.rate_per_kg_day, 'rate_per_kg_day'))


def _lines(values, domain):
    result = tuple(values)
    if any(not isinstance(line, MonoenergeticLine) or line.energy_domain is not domain
           for line in result):
        raise ValueError('all lines must match the declared energy domain')
    return result


def _component_counts(value, bins):
    result = counts(value, 'background response counts')
    if result.shape != (len(bins)-1,):
        raise ValueError('background response counts must match observed bins')
    return result


@dataclass(frozen=True)
class TruthBackground:
    """NR/ER continuum plus lines, folded through this component's response.

    spectrum(E_truth_keV) is events/(kg day keV_nr) or .../keV_er as declared.
    spectrum_energy_range_keV and spectrum_breakpoints_keV apply ONLY to the
    continuum. Every line has its own energy; lines are never quadrature deltas.
    Matrix responses can fold continua but do not provide point-line migration.
    """
    energy_domain: EnergyDomain
    response: TruthResponse
    spectrum: Callable | None = None
    spectrum_energy_range_keV: tuple[float, float] | None = None
    lines: tuple[MonoenergeticLine, ...] = ()
    spectrum_breakpoints_keV: tuple[float, ...] = ()
    provenance: str = 'user-supplied; unverified'

    def __post_init__(self):
        if self.energy_domain not in (EnergyDomain.NR, EnergyDomain.ER):
            raise ValueError('truth background requires explicit EnergyDomain.NR or ER')
        if getattr(self.response, 'input_energy_domain', None) is not self.energy_domain:
            raise ValueError('response input energy domain does not match background')
        lines = _lines(self.lines, self.energy_domain)
        points = array(tuple(self.spectrum_breakpoints_keV), 'spectrum_breakpoints_keV')
        if points.ndim != 1:
            raise ValueError('spectrum breakpoints must be a vector')
        if self.spectrum is None:
            if not lines or self.spectrum_energy_range_keV is not None or points.size:
                raise ValueError('a continuum needs a spectrum/range; a line-only input needs lines')
        else:
            if not callable(self.spectrum) or not callable(getattr(self.response, 'expected_counts', None)):
                raise TypeError('spectrum/response must implement the continuum callable port')
            interval = edges(self.spectrum_energy_range_keV, 'spectrum_energy_range_keV', nonnegative=True)
            if len(interval) != 2:
                raise ValueError('spectrum_energy_range_keV must contain two endpoints')
            object.__setattr__(self, 'spectrum_energy_range_keV', tuple(interval))
        if lines and not callable(getattr(self.response, 'expected_line_counts', None)):
            raise TypeError('response cannot fold point lines; supply a line-capable response')
        object.__setattr__(self, 'lines', lines)
        object.__setattr__(self, 'spectrum_breakpoints_keV', tuple(points))
        object.__setattr__(self, 'provenance', text(self.provenance, 'provenance'))

    def expected_counts(self, observed_bin_edges_keV, exposure_kg_day):
        bins = edges(observed_bin_edges_keV)
        exposure = scalar(exposure_kg_day, 'exposure_kg_day')
        result = np.zeros(len(bins)-1)
        if exposure == 0:
            return result
        with np.errstate(over='raise', invalid='raise'):
            if self.spectrum is not None:
                result += _component_counts(self.response.expected_counts(self.spectrum,
                    bins.copy(), exposure, self.spectrum_energy_range_keV,
                    truth_breakpoints=self.spectrum_breakpoints_keV), bins)
            for line in self.lines:
                result += _component_counts(self.response.expected_line_counts(
                    line.energy_keV, line.rate_per_kg_day, bins.copy(), exposure), bins)
        return output(result)


@dataclass(frozen=True)
class ObservedLineBackground:
    """Already observed point lines; direct bin assignment, no response/efficiency.

    Bins are [a,b), with the final right edge included, as for NoSmearing.
    A measured broadened peak should instead be supplied as an observed
    continuum/histogram; this class does not infer or reapply a line width.
    """
    lines: tuple[MonoenergeticLine, ...]
    provenance: str = 'user-supplied; unverified'

    def __post_init__(self):
        object.__setattr__(self, 'lines', _lines(self.lines, EnergyDomain.OBSERVED))
        object.__setattr__(self, 'provenance', text(self.provenance, 'provenance'))

    def expected_counts(self, observed_bin_edges_keV, exposure_kg_day):
        bins = edges(observed_bin_edges_keV)
        exposure = scalar(exposure_kg_day, 'exposure_kg_day')
        result = np.zeros(len(bins)-1)
        with np.errstate(over='raise', invalid='raise'):
            for line in self.lines:
                if bins[0] <= line.energy_keV <= bins[-1]:
                    j = min(int(np.searchsorted(bins, line.energy_keV, side='right')-1), len(result)-1)
                    result[j] += exposure * line.rate_per_kg_day
        return output(result)
