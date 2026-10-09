"""Shared expected-bin-count validation and explicit integration policy.

No physics models live here. Arrays crossing these interfaces are copied;
configuration containers use tuples. User callbacks must remain pure.
"""
from dataclasses import dataclass
from numbers import Integral, Real
import warnings

import numpy as np
from scipy.integrate import IntegrationWarning, quad


def scalar(value, name, *, positive=False):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError(f'{name} must be a real scalar')
    result = float(value)
    if not np.isfinite(result) or (result <= 0 if positive else result < 0):
        raise ValueError(f'{name} must be finite and {"positive" if positive else "nonnegative"}')
    return result


def array(value, name, *, nonnegative=True):
    result = np.asarray(value)
    if result.dtype.kind not in 'iuf':
        raise ValueError(f'{name} must contain real numbers, excluding booleans')
    result = np.array(result, dtype=float, copy=True)
    if np.any(~np.isfinite(result)) or (nonnegative and np.any(result < 0)):
        raise ValueError(f'{name} has invalid values')
    return result


def edges(value, name='observed_bin_edges_keV', *, nonnegative=False):
    result = array(value, name, nonnegative=nonnegative)
    if result.ndim != 1 or len(result) < 2 or np.any(np.diff(result) <= 0):
        raise ValueError(f'{name} must be a strictly increasing finite vector with at least two edges')
    return result


def counts(value, name='expected counts', *, observed=False):
    result = array(value, name)
    if result.ndim != 1 or result.size == 0:
        raise ValueError(f'{name} must be a nonempty bin-count vector')
    if observed and np.any(result != np.floor(result)):
        raise ValueError('actual observed counts must be nonnegative integers')
    return result


def output(value):
    result = array(value, 'computed output')
    return float(result) if result.ndim == 0 else result


def sample(callback, x, name, *, probability=False):
    value = array(callback(float(x)), name)
    if value.ndim != 0 or (probability and value > 1):
        raise ValueError(f'{name} must return a finite nonnegative scalar'+(' <= 1' if probability else ''))
    return float(value)


def elementwise(callback, x, name, *, probability=False, nonnegative_input=True):
    value = array(x, name+' input', nonnegative=nonnegative_input)
    result = np.empty_like(value)
    for index in np.ndindex(value.shape):
        result[index] = sample(callback, value[index], name, probability=probability)
    return output(result)


def text(value, name):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f'{name} must be a nonempty scalar string')
    return str(value)


@dataclass(frozen=True)
class IntegrationPolicy:
    """Adaptive scalar quadrature policy; epsabs is in the integrated units.

    Known discontinuities/narrow features must be supplied as breakpoints.
    IntegrationWarning is an error; no warning suppression or renormalization.
    """
    epsrel: float = 1e-8
    epsabs: float = 0.0
    limit: int = 300

    def __post_init__(self):
        rel = scalar(self.epsrel, 'epsrel', positive=True)
        if not 50*np.finfo(float).eps < rel < 1:
            raise ValueError('epsrel must lie between 50 machine epsilons and 1')
        if isinstance(self.limit, (bool, np.bool_)) or not isinstance(self.limit, Integral) or self.limit < 1:
            raise ValueError('limit must be a positive integer')
        object.__setattr__(self, 'epsrel', rel)
        object.__setattr__(self, 'epsabs', scalar(self.epsabs, 'epsabs'))
        object.__setattr__(self, 'limit', int(self.limit))

    def integrate(self, callback, lower, upper, *, breakpoints=()):
        bounds = array([lower, upper], 'integration bounds', nonnegative=False)
        if bounds[0] >= bounds[1]:
            raise ValueError('integration bounds must be increasing')
        points = array(tuple(breakpoints), 'breakpoints', nonnegative=False)
        if points.ndim != 1:
            raise ValueError('breakpoints must be a vector')
        points = np.unique(points[(points > lower) & (points < upper)])
        with warnings.catch_warnings(), np.errstate(over='raise', invalid='raise', divide='raise'):
            warnings.simplefilter('error', IntegrationWarning)
            result = quad(callback, float(lower), float(upper), points=points,
                          epsrel=self.epsrel, epsabs=self.epsabs, limit=self.limit)[0]
        return output(result)
