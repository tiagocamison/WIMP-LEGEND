"""Truth differential rates -> expected observed-bin counts, independent of DM."""
from dataclasses import dataclass
from typing import Callable, Protocol
import numpy as np
from scipy.optimize import brentq

from .._analysis import IntegrationPolicy, array, counts, edges, output, sample, scalar, text
from .quenching import ConstantQuenching, LindhardQuenching, QuenchingModel
from .resolution import GaussianResolution, NoSmearing, ResolutionModel
from .efficiency import Efficiency


class CountResponse(Protocol):
    def expected_counts(self, truth_spectrum, observed_bin_edges_keV, exposure_kg_day,
                        truth_energy_range_keV, *, truth_breakpoints=()):
        """Fold events/(kg day keV_nr) to expected observed-energy bin events."""
        ...


def _inputs(truth_spectrum, observed_bin_edges_keV, exposure_kg_day,
            truth_energy_range_keV, truth_breakpoints):
    if not callable(truth_spectrum):
        raise TypeError('truth_spectrum must be callable')
    bins = edges(observed_bin_edges_keV)
    exposure = scalar(exposure_kg_day, 'exposure_kg_day')
    interval = edges(truth_energy_range_keV, 'truth_energy_range_keV', nonnegative=True)
    if len(interval) != 2:
        raise ValueError('truth_energy_range_keV must contain exactly two endpoints')
    points = array(tuple(truth_breakpoints), 'truth_breakpoints')
    if points.ndim != 1:
        raise ValueError('truth_breakpoints must be a vector')
    return bins, exposure, interval, points


def _probabilities(value, n_bins):
    result = array(value, 'migration probabilities')
    if result.shape != (n_bins,) or np.any(result > 1) or np.sum(result) > 1+8*np.finfo(float).eps:
        raise ValueError('migration probabilities must match bins and sum to at most one')
    return result


def _fold(truth_spectrum, bins, exposure, interval, points, kernel, truth_efficiency, policy):
    if not isinstance(policy, IntegrationPolicy):
        raise TypeError('policy must be IntegrationPolicy')
    result = np.zeros(len(bins)-1)
    if exposure == 0:
        return result
    for j in range(len(result)):
        def integrand(energy):
            rate = sample(truth_spectrum, energy, 'truth spectrum')
            efficiency = 1 if truth_efficiency is None else sample(
                truth_efficiency, energy, 'truth efficiency', probability=True)
            probability = _probabilities(kernel(energy), len(result))[j]
            return rate * efficiency * probability
        result[j] = exposure * policy.integrate(integrand, *interval, breakpoints=points)
    return output(result)


@dataclass(frozen=True)
class DetectorResponse:
    quenching: QuenchingModel
    resolution: ResolutionModel
    truth_efficiency: Efficiency | None = None
    observed_efficiency: Efficiency | None = None
    policy: IntegrationPolicy = IntegrationPolicy()
    observed_breakpoints: tuple[float, ...] = ()

    def __post_init__(self):
        if not callable(getattr(self.quenching, 'mean_energy_keV', None)):
            raise TypeError('quenching must implement mean_energy_keV')
        if not callable(getattr(self.resolution, 'bin_probabilities', None)):
            raise TypeError('resolution must implement bin_probabilities')
        for efficiency in (self.truth_efficiency, self.observed_efficiency):
            if efficiency is not None and not callable(efficiency):
                raise TypeError('efficiency must be callable')
        if not isinstance(self.policy, IntegrationPolicy):
            raise TypeError('policy must be IntegrationPolicy')
        points = array(tuple(self.observed_breakpoints), 'observed_breakpoints', nonnegative=False)
        if points.ndim != 1:
            raise ValueError('observed_breakpoints must be a vector')
        object.__setattr__(self, 'observed_breakpoints', tuple(points))

    def bin_probabilities(self, E_nr_keV, observed_bin_edges_keV):
        energy = array(E_nr_keV, 'E_nr_keV')
        mean = array(self.quenching.mean_energy_keV(energy), 'E_ee_keV')
        if mean.shape != energy.shape:
            raise ValueError('quenching must preserve energy shape')
        return self.resolution.bin_probabilities(
            mean, observed_bin_edges_keV, observed_efficiency=self.observed_efficiency,
            policy=self.policy, observed_breakpoints=self.observed_breakpoints)

    def expected_counts(self, truth_spectrum, observed_bin_edges_keV, exposure_kg_day,
                        truth_energy_range_keV, *, truth_breakpoints=()):
        bins, exposure, interval, points = _inputs(truth_spectrum, observed_bin_edges_keV,
            exposure_kg_day, truth_energy_range_keV, truth_breakpoints)
        # Split known migration edges; add width-scale points for a constant
        # Gaussian so quadrature can discover narrow migration windows.
        # These are breakpoints only: the Gaussian is never truncated.
        deterministic = isinstance(self.resolution, NoSmearing) or (
            isinstance(self.resolution, GaussianResolution) and
            not callable(self.resolution.sigma_keV) and self.resolution.sigma_keV == 0)
        constant_gaussian = isinstance(self.resolution, GaussianResolution) and not callable(self.resolution.sigma_keV)
        if (deterministic or constant_gaussian) and isinstance(
                self.quenching, (ConstantQuenching, LindhardQuenching)):
            mean_lo = float(self.quenching.mean_energy_keV(interval[0]))
            mean_hi = float(self.quenching.mean_energy_keV(interval[1]))
            crossings = []
            candidates = bins
            if constant_gaussian and self.resolution.sigma_keV > 0:
                offsets = np.array([-8.,-4.,-2.,-1.,0.,1.,2.,4.,8.])
                candidates = (bins[:,None]+offsets*self.resolution.sigma_keV).ravel()
            for edge in candidates[(candidates > mean_lo) & (candidates < mean_hi)]:
                if isinstance(self.quenching, ConstantQuenching):
                    crossings.append(edge/self.quenching.q)
                else:
                    crossings.append(brentq(lambda energy: self.quenching.mean_energy_keV(energy)-edge,
                                            *interval, xtol=np.nextafter(0.,1.), rtol=1e-14))
            points = np.concatenate((points, crossings))
        return _fold(truth_spectrum, bins, exposure, interval, points,
                     lambda energy: self.bin_probabilities(energy, bins),
                     self.truth_efficiency, self.policy)


@dataclass(frozen=True)
class KernelResponse:
    """User-supplied bin-integrated kernel callback(E_nr_keV, observed_edges).

    includes_efficiency is mandatory. Extra truth efficiency is allowed only
    when false. Observed efficiency must already be integrated by the callback;
    bin-integrated probabilities cannot support a separate sub-bin efficiency.
    """
    kernel: Callable
    includes_efficiency: bool
    truth_efficiency: Efficiency | None = None
    policy: IntegrationPolicy = IntegrationPolicy()

    def __post_init__(self):
        if not callable(self.kernel) or type(self.includes_efficiency) is not bool:
            raise ValueError('kernel must be callable and includes_efficiency an explicit bool')
        if self.includes_efficiency and self.truth_efficiency is not None:
            raise ValueError('kernel already owns efficiencies')
        if self.truth_efficiency is not None and not callable(self.truth_efficiency):
            raise TypeError('truth_efficiency must be callable')
        if not isinstance(self.policy, IntegrationPolicy):
            raise TypeError('policy must be IntegrationPolicy')

    def expected_counts(self, truth_spectrum, observed_bin_edges_keV, exposure_kg_day,
                        truth_energy_range_keV, *, truth_breakpoints=()):
        bins, exposure, interval, points = _inputs(truth_spectrum, observed_bin_edges_keV,
            exposure_kg_day, truth_energy_range_keV, truth_breakpoints)
        return _fold(truth_spectrum, bins, exposure, interval, points,
                     lambda energy: self.kernel(energy, bins.copy()), self.truth_efficiency, self.policy)


@dataclass(frozen=True)
class MatrixResponse:
    """M[j,i] probabilities, observed=M@truth_counts; columns sum <=1.

    Immutable tuple storage. Exact edges and complete truth-bin range required.
    within_bin_averaging is a required description of how M was obtained. This
    piecewise-constant migration model is approximate for varying sub-bin rates.
    """
    truth_bin_edges_keV: tuple[float, ...]
    observed_bin_edges_keV: tuple[float, ...]
    matrix: tuple[tuple[float, ...], ...]
    includes_efficiency: bool
    within_bin_averaging: str
    truth_efficiency: Efficiency | None = None
    policy: IntegrationPolicy = IntegrationPolicy()

    def __post_init__(self):
        truth = edges(self.truth_bin_edges_keV, 'truth_bin_edges_keV', nonnegative=True)
        observed = edges(self.observed_bin_edges_keV)
        matrix = array(self.matrix, 'matrix')
        if matrix.shape != (len(observed)-1, len(truth)-1) or np.any(matrix > 1) or np.any(
                matrix.sum(axis=0) > 1+8*np.finfo(float).eps):
            raise ValueError('matrix shape/probabilities violate M[j,i] contract')
        if type(self.includes_efficiency) is not bool:
            raise ValueError('includes_efficiency must be an explicit bool')
        if self.includes_efficiency and self.truth_efficiency is not None:
            raise ValueError('matrix already owns efficiencies')
        if self.truth_efficiency is not None and not callable(self.truth_efficiency):
            raise TypeError('truth_efficiency must be callable')
        if not isinstance(self.policy, IntegrationPolicy):
            raise TypeError('policy must be IntegrationPolicy')
        object.__setattr__(self, 'truth_bin_edges_keV', tuple(truth))
        object.__setattr__(self, 'observed_bin_edges_keV', tuple(observed))
        object.__setattr__(self, 'matrix', tuple(tuple(row) for row in matrix))
        object.__setattr__(self, 'within_bin_averaging', text(self.within_bin_averaging, 'within_bin_averaging'))

    def apply_truth_counts(self, truth_counts):
        values = counts(truth_counts, 'truth_counts')
        if values.shape != (len(self.truth_bin_edges_keV)-1,):
            raise ValueError('truth counts must match matrix truth bins')
        with np.errstate(over='raise', invalid='raise'):
            return output(np.asarray(self.matrix) @ values)

    def expected_counts(self, truth_spectrum, observed_bin_edges_keV, exposure_kg_day,
                        truth_energy_range_keV, *, truth_breakpoints=()):
        bins, exposure, interval, points = _inputs(truth_spectrum, observed_bin_edges_keV,
            exposure_kg_day, truth_energy_range_keV, truth_breakpoints)
        if not np.array_equal(bins, self.observed_bin_edges_keV) or not np.array_equal(
                interval, (self.truth_bin_edges_keV[0], self.truth_bin_edges_keV[-1])):
            raise ValueError('matrix requires exact observed edges and its complete truth range')
        if exposure == 0:
            return np.zeros(len(bins)-1)
        def integrand(energy):
            efficiency = 1 if self.truth_efficiency is None else sample(
                self.truth_efficiency, energy, 'truth efficiency', probability=True)
            return sample(truth_spectrum, energy, 'truth spectrum') * efficiency
        truth_counts = [exposure*self.policy.integrate(integrand, lo, hi, breakpoints=points)
                        for lo, hi in zip(self.truth_bin_edges_keV[:-1], self.truth_bin_edges_keV[1:])]
        return self.apply_truth_counts(truth_counts)
