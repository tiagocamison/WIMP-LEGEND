"""Conditional reconstructed-energy bin probabilities; no window renormalization."""
from dataclasses import dataclass
from typing import Callable, Protocol
import numpy as np
from scipy.special import ndtr, ndtri
from .._analysis import IntegrationPolicy, array, edges, elementwise, output, sample, scalar
from .efficiency import ConstantEfficiency, Efficiency


class ResolutionModel(Protocol):
    def bin_probabilities(self, E_ee_mean_keV, observed_bin_edges_keV, *,
                          observed_efficiency=None, policy=IntegrationPolicy(),
                          observed_breakpoints=()):
        """Return shape mean.shape+(number_of_bins,), including observed efficiency."""
        ...


def _deterministic(mean, bins, efficiency):
    result = np.zeros(len(bins)-1)
    # Half-open bins [a,b); the last bin includes its right edge.
    if bins[0] <= mean <= bins[-1]:
        j = min(int(np.searchsorted(bins, mean, side='right')-1), len(result)-1)
        result[j] = 1 if efficiency is None else sample(efficiency, mean, 'observed efficiency', probability=True)
    return result


@dataclass(frozen=True)
class NoSmearing:
    def bin_probabilities(self, E_ee_mean_keV, observed_bin_edges_keV, *,
                          observed_efficiency=None, policy=IntegrationPolicy(),
                          observed_breakpoints=()):
        means = array(E_ee_mean_keV, 'E_ee_mean_keV')
        bins = edges(observed_bin_edges_keV)
        result = np.empty(means.shape+(len(bins)-1,))
        for index in np.ndindex(means.shape):
            result[index] = _deterministic(float(means[index]), bins, observed_efficiency)
        return output(result)


@dataclass(frozen=True)
class GaussianResolution:
    """sigma_keV is a constant or scalar callback of quenched E_ee_mean_keV.

    Arbitrary observed efficiency is integrated in probability coordinates
    inside each bin. Positive-tail bins use survival probabilities, avoiding
    cancellation of CDF values close to one. sigma=0 uses NoSmearing semantics.
    """
    sigma_keV: float | Callable

    def __post_init__(self):
        if not callable(self.sigma_keV):
            object.__setattr__(self, 'sigma_keV', scalar(self.sigma_keV, 'sigma_keV'))

    def sigma_at(self, E_ee_mean_keV):
        if callable(self.sigma_keV):
            return elementwise(self.sigma_keV, E_ee_mean_keV, 'sigma_keV')
        energy = array(E_ee_mean_keV, 'E_ee_mean_keV')
        return output(np.full_like(energy, self.sigma_keV))

    def bin_probabilities(self, E_ee_mean_keV, observed_bin_edges_keV, *,
                          observed_efficiency=None, policy=IntegrationPolicy(),
                          observed_breakpoints=()):
        means = array(E_ee_mean_keV, 'E_ee_mean_keV')
        bins = edges(observed_bin_edges_keV)
        points = array(tuple(observed_breakpoints), 'observed_breakpoints', nonnegative=False)
        if points.ndim != 1 or not isinstance(policy, IntegrationPolicy):
            raise ValueError('invalid observed breakpoints or integration policy')
        sigmas = np.asarray(self.sigma_at(means))
        result = np.zeros(means.shape+(len(bins)-1,))
        for index in np.ndindex(means.shape):
            mean, sigma = float(means[index]), float(sigmas[index])
            if sigma == 0:
                result[index] = _deterministic(mean, bins, observed_efficiency)
                continue
            with np.errstate(over='raise', invalid='raise', divide='raise'):
                z = (bins-mean)/sigma
            for j, (a, b) in enumerate(zip(z[:-1], z[1:])):
                positive_tail = a >= 0
                lo, hi = (ndtr(-b), ndtr(-a)) if positive_tail else (ndtr(a), ndtr(b))
                probability = float(hi-lo)
                if observed_efficiency is None or probability == 0:
                    result[index+(j,)] = probability
                elif isinstance(observed_efficiency, ConstantEfficiency):
                    result[index+(j,)] = probability * observed_efficiency.epsilon
                else:
                    def weighted(u):
                        quantile = lo + u*probability
                        coordinate = -ndtri(quantile) if positive_tail else ndtri(quantile)
                        return sample(observed_efficiency, mean+sigma*coordinate,
                                      'observed efficiency', probability=True)
                    p = points[(points > bins[j]) & (points < bins[j+1])]
                    p = ndtr(-(p-mean)/sigma) if positive_tail else ndtr((p-mean)/sigma)
                    p = (p-lo)/probability
                    result[index+(j,)] = probability * policy.integrate(weighted, 0, 1, breakpoints=p)
        return output(result)
