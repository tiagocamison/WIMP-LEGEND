"""Exact known-background Poisson CLs, optionally collapsing bins to totals."""
from dataclasses import dataclass
from math import isfinite, log, log1p
import numpy as np
from scipy.special import gammaln
from scipy.stats import poisson
from .._analysis import scalar
from .base import LimitResult, RootPolicy, confidence, data, total


def _log_tail_sum(n, mean):
    """log(sum_{k=0}^n (n!/k!)*mean^(k-n)), mean>n.

    Positive descending terms are summed until their contribution rounds away.
    This avoids log-CDF underflow and large absolute exponent cancellation.
    """
    term, result = 1., 1.
    for k in range(n, 0, -1):
        term *= k/mean
        previous = result
        result += term
        if result == previous:
            return log(result)
    return log(result)


def _log_cdf(n, mean):
    if mean == 0:
        return 0.
    if n == 0:
        return -mean
    value = float(poisson.logcdf(n, mean))
    if isfinite(value):
        return value
    if mean <= n:
        raise ValueError('Poisson CDF unsupported numerical range')
    return -mean+n*log(mean)-float(gammaln(n+1))+_log_tail_sum(n,mean)


def _log_cls_events(signal_events, background_events, n):
    if signal_events == 0:
        return 0.
    if n == 0:
        # Exact cancellation of known background, even when B+s rounds to B.
        return -signal_events
    mean = background_events + signal_events
    if not isfinite(mean) or mean == background_events:
        raise ValueError('Poisson mean increment is outside supported numerical precision')
    if background_events > n:
        return (-signal_events + n*log1p(signal_events/background_events)
                + _log_tail_sum(n,mean)-_log_tail_sum(n,background_events))
    return _log_cdf(n,mean)-_log_cdf(n,background_events)


@dataclass(frozen=True)
class PoissonCLsCounting:
    """Known-background exact CLs; all supplied bins collapse into one count.

    Default confidence .90 is explicit in each result. Expected sensitivity is
    the smallest-integer Poisson .5 quantile, not n=B or a rounded Asimov count.
    No nuisance uncertainty or bin shape information is used.
    """
    root_policy: RootPolicy = RootPolicy()

    def __post_init__(self):
        if not isinstance(self.root_policy, RootPolicy):
            raise TypeError('root_policy must be RootPolicy')

    def log_cls(self, mu, signal, background, observed):
        s,b,n = data(signal,background,observed)
        strength = scalar(mu, 'mu')
        signal_total = total(s)
        events = strength*signal_total
        if not isfinite(events) or (strength > 0 and signal_total > 0 and events == 0):
            raise ValueError('signal event total is outside floating-point range')
        n_total = total(n)
        if n_total >= 2**53:
            raise ValueError('total observed count exceeds exact integer precision')
        return _log_cls_events(events,total(b),int(n_total))

    def acceptance_margin(self, signal, background, observed, *, confidence_level=.90):
        cl = confidence(confidence_level)
        return self.log_cls(1.,signal,background,observed)-log1p(-cl)

    def acceptance_metadata(self):
        return ('poisson_cls_counting', 'log CLs(mu=1)-log(1-CL); known B; bins collapsed',
                'exact discrete Poisson construction')

    def _limit(self, signal, background, observed, cl, kind, notes=()):
        s,b,n = data(signal,background,observed)
        signal_total, bg, n_total = total(s),total(b),total(n)
        if n_total >= 2**53:
            raise ValueError('total observed count exceeds exact integer precision')
        target = log1p(-cl)
        if signal_total == 0:
            limit, status = np.inf, 'no_sensitivity'
        else:
            # Solve in signal events, then convert to strength. This avoids
            # bracket dependence on template normalization.
            if n_total == 0:
                events = -target
            else:
                events = self.root_policy.crossing(
                    lambda x: target-_log_cls_events(x,bg,int(n_total)),
                    initial=max(1.,n_total**.5))
            limit = events/signal_total
            if not isfinite(limit) or limit == 0:
                raise ValueError('signal strength limit is outside numerical range')
            status = 'success'
        return LimitResult(limit,'poisson_cls_counting',
            'P(N<=n|B+mu*S) / P(N<=n|B); known B; bins collapsed',cl,kind,
            'exact discrete Poisson construction; numerical root',status,True,
            tuple(n),tuple(notes))

    def observed_upper_limit(self, signal, background, observed, *, confidence_level=.90):
        return self._limit(signal,background,observed,confidence(confidence_level),'observed')

    def expected_upper_limit(self, signal, background, *, confidence_level=.90):
        s,b,_ = data(signal,background)
        bg = total(b)
        if bg >= 2**52:
            raise ValueError('expected Poisson quantile outside supported integer precision')
        n = float(poisson.ppf(.5,bg))
        if not isfinite(n):
            raise ValueError('Poisson median quantile failed')
        # Check the declared discrete quantile, including means near a jump.
        if poisson.cdf(n,bg) < .5 or (n > 0 and poisson.cdf(n-1,bg) >= .5):
            raise ValueError('Poisson median quantile failed its discrete bracket')
        return self._limit([total(s)],[bg],[n],confidence(confidence_level),'expected_median',
                           ('N~Poisson(B); smallest integer n with CDF(n)>=0.5',))
