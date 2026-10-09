"""Replaceable count-only statistical contracts and auditable numerical results."""
from dataclasses import dataclass
from math import fsum, isfinite
from numbers import Integral
from typing import Protocol
import numpy as np
from scipy.optimize import brentq
from .._analysis import counts, scalar, text


def data(signal, background, observed=None, *, asimov=False):
    s, b = counts(signal, 'signal template'), counts(background, 'background template')
    if s.shape != b.shape:
        raise ValueError('signal and background must have identical shapes')
    n = None if observed is None else counts(observed, 'observed counts', observed=not asimov)
    if n is not None and n.shape != s.shape:
        raise ValueError('observed counts must match expected-count shapes')
    if n is not None and not asimov and np.any(n >= 2**53):
        raise ValueError('observed counts exceed exact integer float precision')
    return s, b, n


def total(values):
    try:
        result = fsum(float(value) for value in values)
    except OverflowError as exc:
        raise ValueError('total expected counts overflow') from exc
    if not isfinite(result):
        raise ValueError('total counts are not finite')
    return result


def confidence(value, *, asymptotic=False):
    value = scalar(value, 'confidence_level', positive=True)
    if not (0.5 if asymptotic else 0) < value < 1:
        raise ValueError('confidence_level outside supported interval')
    return value


@dataclass(frozen=True)
class RootPolicy:
    xtol: float = 1e-10
    rtol: float = 1e-10
    maxiter: int = 200
    max_bracket_steps: int = 128

    def __post_init__(self):
        for name in ('xtol', 'rtol'):
            object.__setattr__(self, name, scalar(getattr(self, name), name, positive=True))
        if self.rtol < 4*np.finfo(float).eps:
            raise ValueError('rtol must be at least four machine epsilons')
        for name in ('maxiter', 'max_bracket_steps'):
            value = getattr(self, name)
            if isinstance(value, (bool, np.bool_)) or not isinstance(value, Integral) or value < 1:
                raise ValueError(f'{name} must be a positive integer')
            object.__setattr__(self, name, int(value))

    def crossing(self, function, *, lower=0., initial=1.):
        """Root of a bracketed increasing sign change; failures are explicit."""
        lower = scalar(lower, 'root lower bound')
        initial = scalar(initial, 'root initial step', positive=True)
        left = float(function(lower))
        if np.isnan(left) or left > 0:
            raise ValueError('root lower bound does not bracket the accepted region')
        if left == 0:
            return lower
        upper = lower + initial
        for _ in range(self.max_bracket_steps):
            if not isfinite(upper) or upper <= lower:
                raise ValueError('root bracket is outside numerical range')
            right = float(function(upper))
            if np.isnan(right):
                raise ValueError('root function returned NaN')
            if right >= 0:
                break
            upper = lower + 2*(upper-lower)
        else:
            raise RuntimeError('failed to bracket statistical root')
        root, result = brentq(function, lower, upper, xtol=self.xtol,
                              rtol=self.rtol, maxiter=self.maxiter, full_output=True, disp=False)
        if not result.converged or not isfinite(root):
            raise RuntimeError('statistical root did not converge')
        return float(root)


@dataclass(frozen=True)
class LimitResult:
    mu_upper: float
    method: str
    construction: str
    confidence_level: float
    kind: str
    approximation: str
    status: str
    converged: bool
    reference_counts: tuple[float, ...] = ()
    notes: tuple[str, ...] = ()
    strength_convention: str = 'lambda_i = mu * S_i + B_i; S_i supplied at mu=1'

    def __post_init__(self):
        # Infinity is reserved for an explicitly insensitive zero template.
        if self.mu_upper == np.inf and self.status != 'no_sensitivity':
            raise ValueError('infinite limit requires no_sensitivity status')
        if self.mu_upper != np.inf:
            object.__setattr__(self, 'mu_upper', scalar(self.mu_upper, 'mu_upper'))
        object.__setattr__(self, 'confidence_level', confidence(self.confidence_level))
        for name in ('method','construction','kind','approximation','status','strength_convention'):
            object.__setattr__(self, name, text(getattr(self,name), name))
        if type(self.converged) is not bool:
            raise ValueError('converged must be bool')
        reference = tuple(scalar(x, 'reference count') for x in self.reference_counts)
        object.__setattr__(self, 'reference_counts', reference)
        object.__setattr__(self, 'notes', tuple(text(x, 'note') for x in self.notes))


class StatisticsBackend(Protocol):
    def observed_upper_limit(self, signal, background, observed, *, confidence_level=.90) -> LimitResult:
        ...

    def expected_upper_limit(self, signal, background, *, confidence_level=.90) -> LimitResult:
        ...

    def acceptance_margin(self, signal, background, observed, *, confidence_level=.90) -> float:
        """Nonnegative means accepted; test supplied absolute signal at strength 1."""
        ...

    def acceptance_metadata(self) -> tuple[str, str, str]:
        """Return method ID, test construction, and approximation status for scans."""
        ...
