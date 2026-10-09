"""Independent-bin Poisson shape likelihood and explicitly asymptotic limits."""
from dataclasses import dataclass
from math import isfinite, log1p
import numpy as np
from scipy.special import gammaln, xlogy
from scipy.stats import norm
from .._analysis import scalar
from .base import LimitResult, RootPolicy, confidence, data, total


def _means(mu,s,b):
    mu = scalar(mu,'mu')
    with np.errstate(over='raise', invalid='raise'):
        means = mu*s+b
    if np.any(~np.isfinite(means)) or (mu > 0 and np.any((s > 0) & (means == b))):
        raise ValueError('signal mean increment is outside supported numerical precision')
    return means


def _log_likelihood(mu,s,b,n):
    means = _means(mu,s,b)
    if np.any((means == 0) & (n > 0)):
        return -np.inf
    with np.errstate(over='raise', invalid='raise'):
        value = np.sum(xlogy(n,means)-means-gammaln(n+1))
    if not isfinite(value):
        raise ValueError('log likelihood outside numerical range')
    return float(value)


@dataclass(frozen=True)
class FitResult:
    mu_hat: float
    log_likelihood: float
    converged: bool
    status: str


@dataclass(frozen=True)
class BinnedPoissonProfile:
    """Known-background one-sided profile limits with q_mu=z_CL^2.

    Asymptotic only; no exact coverage or sparse-count quantitative calibration.
    The separate Asimov path uses fractional n=B explicitly. Actual observations
    must be integer counts. No nuisance profiling or toy Monte Carlo is supplied.
    """
    root_policy: RootPolicy = RootPolicy()

    def __post_init__(self):
        if not isinstance(self.root_policy,RootPolicy):
            raise TypeError('root_policy must be RootPolicy')

    def log_likelihood(self, mu, signal, background, observed):
        s,b,n = data(signal,background,observed)
        return _log_likelihood(mu,s,b,n)

    def _fit(self,s,b,n):
        if np.any((n > 0) & (s == 0) & (b == 0)):
            raise ValueError('observations have zero probability for every signal strength')
        strength = total(s)
        if strength == 0:
            return FitResult(0.,_log_likelihood(0,s,b,n),True,'no_sensitivity')
        def derivative(mu):
            means = _means(mu,s,b)
            if np.any((means == 0) & (s > 0) & (n > 0)):
                return np.inf
            positive = means > 0
            with np.errstate(over='raise', invalid='raise', divide='raise'):
                terms = np.zeros_like(s)
                terms[positive] = s[positive]*n[positive]/means[positive]
            return total(terms)-strength
        if derivative(0.) <= 0:
            mu = 0.
        else:
            mu = self.root_policy.crossing(lambda x:-derivative(x),
                                           initial=max(1.,total(n))/strength)
        return FitResult(mu,_log_likelihood(mu,s,b,n),True,'success')

    def fit(self, signal, background, observed):
        return self._fit(*data(signal,background,observed))

    def _q(self,mu,s,b,n,fit):
        mu = scalar(mu,'mu')
        if fit.mu_hat > mu:
            return 0.
        if mu == fit.mu_hat:
            return 0.
        means, best = _means(mu,s,b),_means(fit.mu_hat,s,b)
        terms = []
        # Direct likelihood ratio with log1p; no subtraction of large log Ls.
        for si,ni,mi,bi in zip(s,n,means,best):
            delta = (mu-fit.mu_hat)*si
            if bi == 0:
                if ni > 0:
                    raise ValueError('profile best fit has zero likelihood')
                terms.append(float(mi))
            elif mi == 0 and ni > 0:
                return np.inf
            else:
                terms.append(float(delta-ni*log1p(delta/bi)))
        from math import fsum
        value = 2*fsum(terms)
        if np.isnan(value) or value < 0:
            raise ValueError('profile likelihood ratio lost numerical precision')
        return value

    def q_mu(self, mu, signal, background, observed):
        s,b,n = data(signal,background,observed)
        return self._q(mu,s,b,n,self._fit(s,b,n))

    def acceptance_margin(self, signal, background, observed, *, confidence_level=.90):
        cl = confidence(confidence_level,asymptotic=True)
        return float(norm.ppf(cl))**2-self.q_mu(1.,signal,background,observed)

    def acceptance_metadata(self):
        return ('binned_poisson_profile', 'z_CL^2-q_mu(mu=1); alternative fit along supplied shape',
                'asymptotic; sparse-count quantitative coverage unsupported/not validated')

    def _limit(self,s,b,n,cl,kind):
        fit = self._fit(s,b,n)
        strength = total(s)
        if strength == 0:
            limit, status = np.inf,'no_sensitivity'
        else:
            threshold = float(norm.ppf(cl))**2
            limit = self.root_policy.crossing(lambda x:self._q(x,s,b,n,fit)-threshold,
                lower=fit.mu_hat,initial=max(1.,total(n)**.5)/strength)
            status = 'asymptotic_unvalidated'
        return LimitResult(limit,'binned_poisson_profile',
            'constrained mu_hat>=0; q_mu=0 when mu_hat>mu; q_mu=z_CL^2',cl,kind,
            'asymptotic one-sided profile likelihood',status,True,tuple(n),
            ('Sparse-count quantitative coverage unsupported/not validated; no nuisances or toys',))

    def observed_upper_limit(self, signal, background, observed, *, confidence_level=.90):
        s,b,n = data(signal,background,observed)
        return self._limit(s,b,n,confidence(confidence_level,asymptotic=True),'observed')

    def expected_upper_limit(self, signal, background, *, confidence_level=.90):
        s,b,_ = data(signal,background)
        return self._limit(s,b,b.copy(),confidence(confidence_level,asymptotic=True),'expected_asimov')
