"""Callback-based binned Poisson profiling, independent of the physics model."""
from dataclasses import dataclass
from numbers import Integral, Real
from types import MappingProxyType
from typing import Mapping

import numpy as np
from scipy.linalg import solve_triangular
from scipy.optimize import minimize
from scipy.stats import norm

from .._analysis import array, counts, scalar, text
from .base import LimitResult, RootPolicy, confidence


def _real(value, name, *, allow_inf=False):
    if (isinstance(value, (bool, np.bool_)) or not isinstance(value, Real)
            or not (np.isfinite(value) or allow_inf and value == np.inf)):
        raise ValueError(f'{name} must be a finite real scalar'
                         + (' or +infinity' if allow_inf else ''))
    return float(value)


@dataclass(frozen=True)
class NuisanceParameter:
    name: str
    nominal: float = 0.
    initial: float | None = None
    bounds: tuple[float | None, float | None] = (None, None)

    def __post_init__(self):
        name = text(self.name, 'parameter name')
        nominal = _real(self.nominal, 'nominal value')
        initial = nominal if self.initial is None else _real(self.initial, 'initial value')
        lo, hi = (None if x is None else _real(x, 'bound') for x in self.bounds)
        if lo is not None and hi is not None and lo > hi:
            raise ValueError('bounds must be ordered')
        for x in (nominal, initial):
            if (lo is not None and x < lo) or (hi is not None and x > hi):
                raise ValueError(f'{name}: nominal and initial values must lie within bounds')
        for key, value in dict(name=name, nominal=nominal, initial=initial, bounds=(lo, hi)).items():
            object.__setattr__(self, key, value)


@dataclass(frozen=True, init=False)
class GaussianPenalty:
    """Mahalanobis squared distance in -2 log L units; covariance follows center order."""
    names: tuple[str, ...]
    center: tuple[float, ...]
    _factor: tuple[tuple[float, ...], ...]

    def __init__(self, center: Mapping[str, float], covariance):
        names = tuple(text(k, 'Gaussian parameter name') for k in center)
        means = tuple(_real(v, 'Gaussian center') for v in center.values())
        cov = array(covariance, 'covariance', nonnegative=False)
        if (not names or cov.shape != (len(names), len(names))
                or not np.array_equal(cov, cov.T)):
            raise ValueError('covariance must be symmetric and match the nonempty center')
        try:
            factor = np.linalg.cholesky(cov)
        except np.linalg.LinAlgError as exc:
            raise ValueError('covariance must be positive definite') from exc
        object.__setattr__(self, 'names', names)
        object.__setattr__(self, 'center', means)
        object.__setattr__(self, '_factor', tuple(map(tuple, factor)))

    @classmethod
    def independent(cls, center, sigma):
        if set(center) != set(sigma):
            raise ValueError('sigma must cover exactly the Gaussian center names')
        widths = np.array([scalar(sigma[k], 'Gaussian sigma', positive=True) for k in center])
        return cls(center, np.diag(widths**2))

    def __call__(self, theta):
        residual = np.array([theta[k] for k in self.names])-self.center
        whitened = solve_triangular(self._factor, residual, lower=True)
        return float(whitened @ whitened)


@dataclass(frozen=True)
class ProfileFitResult:
    mu_hat: float
    theta_hat: Mapping[str, float]
    deviance: float
    converged: bool
    status: str
    evaluations: int
    locally_flat: tuple[str, ...] = ()


def _poisson_deviance(n, means):
    """Twice the log likelihood relative to the saturated Poisson model."""
    if np.any((n > 0) & (means == 0)):
        return np.inf
    positive = n > 0
    y, m = n[positive], means[positive]
    delta = m-y
    near = abs(delta) < .01*y
    terms = np.empty_like(y)
    x = delta[near]/y[near]
    # x-log1p(x), with cancellation avoided close to the minimum.
    tiny = abs(x) < 1e-4
    local = np.empty_like(x)
    t = x[tiny]
    local[tiny] = t*t*(.5+t*(-1/3+t*(.25-t/5)))
    local[~tiny] = x[~tiny]-np.log1p(x[~tiny])
    terms[near] = y[near]*local
    terms[~near] = delta[~near]+y[~near]*(np.log(y[~near])-np.log(m[~near]))
    return float(2*(np.sum(terms)+np.sum(means[~positive])))


@dataclass(frozen=True)
class ProfiledBinnedPoisson:
    """Fit expected_counts(mu, theta) over all supplied bins.

    Penalty is in -2 log L units. +inf marks excluded points. Expected
    counts must be finite and nonnegative (None may mark an invalid point).
    The optimizer is local; supply extra starts for multimodal models.
    """
    expected_counts: object
    parameters: tuple[NuisanceParameter, ...] = ()
    penalty: object = None
    root_policy: RootPolicy = RootPolicy()
    mu_initial: float = 1.
    maxiter: int = 2000

    def __post_init__(self):
        if not callable(self.expected_counts) or (self.penalty is not None and not callable(self.penalty)):
            raise TypeError('expected_counts and any penalty must be callable')
        params = tuple(self.parameters)
        if not all(isinstance(p, NuisanceParameter) for p in params):
            raise TypeError('parameters must contain NuisanceParameter objects')
        names = {p.name for p in params}
        if len(names) != len(params):
            raise ValueError('nuisance parameter names must be unique')
        if isinstance(self.penalty, GaussianPenalty) and not set(self.penalty.names) <= names:
            raise ValueError('Gaussian penalty names must be declared nuisance parameters')
        if not isinstance(self.root_policy, RootPolicy):
            raise TypeError('root_policy must be RootPolicy')
        if isinstance(self.maxiter, (bool, np.bool_)) or not isinstance(self.maxiter, Integral) or self.maxiter < 1:
            raise ValueError('maxiter must be a positive integer')
        object.__setattr__(self, 'parameters', params)
        object.__setattr__(self, 'mu_initial', scalar(self.mu_initial, 'mu_initial', positive=True))

    def _theta(self, values):
        return MappingProxyType(dict(zip((p.name for p in self.parameters), map(float, values))))

    def _means(self, mu, theta):
        result = self.expected_counts(mu, theta)
        return None if result is None else counts(result, 'expected counts')

    def _objective(self, n, mu, theta):
        means = self._means(mu, theta)
        if means is None:
            return np.inf
        if means.shape != n.shape:
            raise ValueError('observed and expected bin counts must have identical shapes')
        penalty = 0. if self.penalty is None else _real(self.penalty(theta), 'penalty', allow_inf=True)
        return _poisson_deviance(n, means)+penalty

    @staticmethod
    def _observed(observed):
        n = counts(observed, 'observed counts', observed=True)
        if np.any(n >= 2**53):
            raise ValueError('observed counts exceed exact integer float precision')
        return n

    def _fit(self, n, mu, starts, warm=None):
        free = tuple(i for i,p in enumerate(self.parameters)
                     if p.bounds[0] is None or p.bounds[0] != p.bounds[1])
        bounds = ([(0., None)] if mu is None else []) + [self.parameters[i].bounds for i in free]
        initial = [p.initial for p in self.parameters]
        nominal = [p.nominal for p in self.parameters]
        seeds = [(self.mu_initial, initial), (0., nominal), (self.mu_initial, nominal)]
        if warm is not None:
            seeds.append((warm.mu_hat, tuple(warm.theta_hat.values())))
        for strength, theta in starts:
            if set(theta) != {p.name for p in self.parameters}:
                raise ValueError('extra start must specify every nuisance parameter')
            values = [_real(theta[p.name], p.name) for p in self.parameters]
            for value,p in zip(values,self.parameters):
                lo,hi = p.bounds
                if (lo is not None and value < lo) or (hi is not None and value > hi):
                    raise ValueError(f'{p.name}: start outside bounds')
            seeds.append((scalar(strength, 'start mu'), values))

        def unpack(x):
            values = nominal.copy()
            for i,value in zip(free, x if mu is not None else x[1:]):
                values[i] = value
            return float(mu if mu is not None else x[0]), self._theta(values)

        evaluations = 0
        def objective(x):
            nonlocal evaluations
            evaluations += 1
            return self._objective(n, *unpack(x))

        if not bounds:
            value = objective(())
            return ProfileFitResult(*unpack(()), value, bool(np.isfinite(value)),
                                    'fixed' if np.isfinite(value) else 'no finite likelihood', evaluations)
        sentinel = 1e100

        def finite_objective(x):
            value = objective(x)
            if np.isposinf(value):
                return sentinel
            if not np.isfinite(value) or abs(value) >= sentinel:
                raise ValueError('objective outside supported numerical range')
            return value

        best, value, success, message = None, np.inf, False, 'no finite likelihood'
        visited = set()
        for strength, theta in seeds:
            x = tuple(([strength] if mu is None else []) + [theta[i] for i in free])
            if x in visited:
                continue
            visited.add(x)
            seed_value = finite_objective(x)
            result = minimize(finite_objective, x, method='Powell', bounds=bounds,
                              options={'maxiter': self.maxiter, 'xtol': 1e-10, 'ftol': 1e-12})
            candidate, loss = (result.x, result.fun) if result.fun < seed_value else (x, seed_value)
            ok = bool(result.success and loss < sentinel and result.fun-loss <= 1e-8)
            if best is None or loss < value or (loss == value and ok and not success):
                best, value, success = np.asarray(candidate), float(loss), ok
                message = str(result.message) if loss < sentinel else 'no finite likelihood'
        flat = []
        labels = (['mu'] if mu is None else []) + [self.parameters[i].name for i in free]
        if success:
            for j, (lo, hi) in enumerate(bounds):
                probes = []
                for direction in (-1., 1.):
                    point = best.copy()
                    step = .01*max(1., abs(point[j]))
                    point[j] = np.clip(point[j]+direction*step,
                                       -np.inf if lo is None else lo, np.inf if hi is None else hi)
                    if point[j] != best[j]:
                        probes.append(objective(point))
                if probes and all(abs(v-value) <= 1e-10 for v in probes):
                    flat.append(labels[j])
        return ProfileFitResult(*unpack(best), value if value < sentinel else np.inf,
                                success, message, evaluations, tuple(flat))

    def fit(self, observed, *, mu=None, starts=()):
        """Unconditional if mu is omitted, otherwise conditional at fixed mu.

        Extra starts are (mu, mapping) pairs; the start mu is ignored in
        conditional fits. Flat-coordinate diagnostics do not detect all ridges.
        """
        fixed = None if mu is None else scalar(mu, 'mu')
        return self._fit(self._observed(observed), fixed, starts)

    @staticmethod
    def _require(fit):
        if not fit.converged:
            raise RuntimeError(f'profile fit failed: {fit.status}')

    def _q(self, mu, n, best, starts):
        self._require(best)
        if mu <= best.mu_hat:
            return 0.
        conditional = self._fit(n, mu, starts, warm=best)
        self._require(conditional)
        q = conditional.deviance-best.deviance
        if q < -1e-8:
            raise RuntimeError('conditional fit improved on global fit; provide better starts')
        return max(0., q)

    def q_mu(self, mu, observed, *, starts=()):
        mu = scalar(mu, 'mu')
        n, starts = self._observed(observed), tuple(starts)
        return self._q(mu, n, self._fit(n, None, starts), starts)

    def _limit(self, n, cl, kind, starts):
        cl, starts = confidence(cl, asymptotic=True), tuple(starts)
        best = self._fit(n, None, starts)
        self._require(best)
        threshold = float(norm.ppf(cl))**2
        limit = self.root_policy.crossing(lambda mu: self._q(mu, n, best, starts)-threshold,
                                          lower=best.mu_hat, initial=self.mu_initial)
        return LimitResult(limit, 'binned_poisson_profile_nuisance',
                           'one-sided q_mu=z_CL^2, mu_hat>=0', cl, kind,
                           'asymptotic one-sided profile likelihood', 'asymptotic_unvalidated',
                           True, tuple(n),
                           ('Not CL_s; requires a regular model and a single upper crossing.',),
                           'user-supplied expected_counts(mu, theta)')

    def observed_upper_limit(self, observed, *, confidence_level=.90, starts=()):
        return self._limit(self._observed(observed), confidence_level, 'observed', starts)

    def expected_upper_limit(self, *, confidence_level=.90, starts=()):
        """Use the fractional Asimov counts at mu=0 and nominal nuisance values."""
        n = self._means(0., self._theta(p.nominal for p in self.parameters))
        if n is None:
            raise ValueError('nominal background-only model is outside the count domain')
        return self._limit(n, confidence_level, 'expected_asimov', starts)
