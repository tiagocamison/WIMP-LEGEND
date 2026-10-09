"""Independent Poisson sums, elementary likelihoods, and synthetic shape tests."""
from dataclasses import FrozenInstanceError
from math import exp, factorial, lgamma, log, sqrt
import numpy as np
import pytest
from scipy.optimize import brentq
from scipy.stats import norm
from wimp_legend.statistics import PoissonCLsCounting, BinnedPoissonProfile, RootPolicy


def poisson_cdf_sum(n,mean):
    return exp(-mean)*sum(mean**k/factorial(k) for k in range(n+1))


@pytest.mark.parametrize('background',[0.,.3,12.,1000.,1e200])
def test_exact_zero_observation_limit_known_background(background):
    result=PoissonCLsCounting().observed_upper_limit([1],[background],[0])
    assert result.mu_upper==pytest.approx(-log(.1),rel=2e-15)
    assert result.method=='poisson_cls_counting'
    assert result.confidence_level==.9 and result.kind=='observed' and result.converged
    assert 'exact' in result.approximation
    assert PoissonCLsCounting().log_cls(2,[1],[background],[0])==-2


@pytest.mark.parametrize('n,bg',[(1,0.),(4,1.3),(6,5.),(20,25.)])
def test_exact_cls_nonzero_observations_independent_sum(n,bg):
    backend=PoissonCLsCounting()
    root=brentq(lambda s:poisson_cdf_sum(n,bg+s)/poisson_cdf_sum(n,bg)-.1,0,100)
    result=backend.observed_upper_limit([2],[bg],[n])
    assert result.mu_upper==pytest.approx(root/2,rel=2e-10)
    values=[backend.log_cls(x,[2],[bg],[n]) for x in [0,.1,1,5,10]]
    assert values[0]==0 and np.all(np.diff(values)<0)
    assert backend.acceptance_margin([2*result.mu_upper],[bg],[n])==pytest.approx(0,abs=2e-9)


@pytest.mark.parametrize('background',[0.,.001,.5,.8,3.2,15.])
def test_expected_median_discrete_poisson(background):
    # Find smallest count with cumulative mass >= 1/2 from elementary probabilities.
    cumulative=0.;n=0
    while True:
        cumulative+=exp(-background)*background**n/factorial(n)
        if cumulative>=.5:break
        n+=1
    result=PoissonCLsCounting().expected_upper_limit([1],[background])
    assert result.reference_counts==(float(n),)
    root=brentq(lambda s:poisson_cdf_sum(n,background+s)/poisson_cdf_sum(n,background)-.1,0,100)
    assert result.mu_upper==pytest.approx(root,rel=2e-10)
    assert result.kind=='expected_median'
    assert 'smallest integer' in result.notes[0]


def test_large_background_logcdf_and_template_scale():
    backend=PoissonCLsCounting()
    # n=1: ratio=e^-s*(1+B+s)/(1+B), even when each CDF underflows.
    bg=1e5
    root=brentq(lambda s:-s+log((1+bg+s)/(1+bg))-log(.1),0,10)
    assert backend.observed_upper_limit([1],[bg],[1]).mu_upper==pytest.approx(root,rel=2e-10)
    high=backend.expected_upper_limit([1],[1e6])
    assert high.mu_upper>1000 and np.isfinite(high.mu_upper)
    tiny=backend.observed_upper_limit([1e-150],[0],[0])
    assert tiny.mu_upper==pytest.approx(-log(.1)*1e150)
    assert backend.observed_upper_limit([1,2],[3,4],[1,2]).mu_upper==pytest.approx(
        backend.observed_upper_limit([3],[7],[3]).mu_upper)


@pytest.mark.parametrize('backend',[PoissonCLsCounting(),BinnedPoissonProfile()])
def test_zero_signal_explicit_no_sensitivity(backend):
    result=backend.observed_upper_limit([0,0],[0,2],[0,1])
    assert result.mu_upper==np.inf and result.status=='no_sensitivity' and result.converged
    assert backend.expected_upper_limit([0],[0]).mu_upper==np.inf
    with pytest.raises(FrozenInstanceError):result.method='other'


@pytest.mark.parametrize('n,mean',[(0,0),(0,3),(2,3),(9,7)])
def test_one_bin_log_likelihood(n,mean):
    got=BinnedPoissonProfile().log_likelihood(0,[1],[mean],[n])
    expected=0 if mean==0 and n==0 else n*log(mean)-mean-lgamma(n+1)
    assert got==pytest.approx(expected,abs=2e-14)
    assert BinnedPoissonProfile().log_likelihood(0,[1],[0],[1])==-np.inf


@pytest.mark.parametrize('n,bg,signal',[(12,3.,2.),(2,3.,2.),(0,0.,1.),(9,0.,3.)])
def test_one_bin_constrained_fit_and_likelihood_ratio(n,bg,signal):
    backend=BinnedPoissonProfile()
    muhat=max(0,(n-bg)/signal)
    fit=backend.fit([signal],[bg],[n])
    assert fit.mu_hat==pytest.approx(muhat,abs=2e-10) and fit.converged
    for mu in [0.,.2,2.,10.]:
        if mu<muhat:expected=0.
        elif mu==muhat:expected=0.
        else:
            best=bg+muhat*signal;tested=bg+mu*signal
            expected=2*(tested-best-(n*log(tested/best) if n else 0))
        assert backend.q_mu(mu,[signal],[bg],[n])==pytest.approx(expected,rel=2e-10,abs=2e-11)


def test_shape_independent_ratio_fit_and_one_sided_threshold():
    backend=BinnedPoissonProfile()
    signal=np.array([10.,40.]);background=np.array([100.,200.]);n=np.array([130,230])
    fit=backend.fit(signal,background,n)
    oracle=brentq(lambda m:sum(s*(ni/(m*s+b)-1) for s,b,ni in zip(signal,background,n)),0,10)
    assert fit.mu_hat==pytest.approx(oracle,rel=1e-10)
    mu=2.
    expected=2*sum((mu-oracle)*s-ni*log((mu*s+b)/(oracle*s+b)) for s,b,ni in zip(signal,background,n))
    assert backend.q_mu(mu,signal,background,n)==pytest.approx(expected,rel=2e-10)
    result=backend.observed_upper_limit(signal,background,n)
    assert backend.q_mu(result.mu_upper,signal,background,n)==pytest.approx(norm.ppf(.9)**2,rel=2e-10)
    assert result.status=='asymptotic_unvalidated' and 'Sparse-count' in result.notes[0]


def test_identical_shape_splitting_no_spurious_information():
    backend=BinnedPoissonProfile()
    collapsed=backend.observed_upper_limit([30],[300],[330])
    split=backend.observed_upper_limit([10,20],[100,200],[110,220])
    assert split.mu_upper==pytest.approx(collapsed.mu_upper,rel=1e-10)
    assert backend.q_mu(2,[10,20],[100,200],[110,220])==pytest.approx(
        backend.q_mu(2,[30],[300],[330]),rel=1e-11)


def test_high_count_shape_improves_and_distinct_templates():
    backend=BinnedPoissonProfile()
    s=[100,0];b=[1000,9000];n=[1000,9000]
    shape=backend.observed_upper_limit(s,b,n)
    collapsed=backend.observed_upper_limit([100],[10000],[10000])
    assert shape.mu_upper<collapsed.mu_upper/2.5
    expected=brentq(lambda m:2*(m*100-1000*log(1+m*.1))-norm.ppf(.9)**2,1e-5,10)
    assert shape.mu_upper==pytest.approx(expected,rel=1e-10)
    assert backend.log_likelihood(1,[100,0],b,n)!=backend.log_likelihood(1,[0,100],b,n)
    assert shape.mu_upper<PoissonCLsCounting().observed_upper_limit(s,b,n).mu_upper


def test_asimov_fractional_counts_explicit_path():
    backend=BinnedPoissonProfile()
    s=[1.,2.];b=[.2,1.7]
    result=backend.expected_upper_limit(s,b)
    expected=brentq(lambda m:2*sum(m*si-bi*log(1+m*si/bi) for si,bi in zip(s,b))-
                    norm.ppf(.9)**2,1e-6,20)
    assert result.mu_upper==pytest.approx(expected,rel=1e-10)
    assert result.kind=='expected_asimov' and result.reference_counts==tuple(b)
    with pytest.raises(ValueError):backend.log_likelihood(1,s,b,b)
    with pytest.raises(ValueError):backend.fit([0],[0],[1])


@pytest.mark.parametrize('backend',[PoissonCLsCounting(),BinnedPoissonProfile()])
@pytest.mark.parametrize('s,b,n',[
    ([-1],[0],[0]),([np.nan],[0],[0]),([1],[np.inf],[0]),([1],[0],[-1]),
    ([1],[0],[.2]),([1],[0],[True]),([1,2],[0],[0]),([1],[0],[0,1]),
    ([],[],[]),([True],[0],[0]),([1],[0],[2**53]),
])
def test_invalid_statistical_data(backend,s,b,n):
    with pytest.raises(ValueError):backend.observed_upper_limit(s,b,n)


@pytest.mark.parametrize('cl',[0,1,-1,np.nan,np.inf,True])
def test_invalid_confidence(cl):
    for backend in [PoissonCLsCounting(),BinnedPoissonProfile()]:
        with pytest.raises(ValueError):backend.expected_upper_limit([1],[1],confidence_level=cl)
    with pytest.raises(ValueError):BinnedPoissonProfile().expected_upper_limit([1],[1],confidence_level=.5)


def test_root_failures_and_data_ownership():
    with pytest.raises(RuntimeError):RootPolicy(max_bracket_steps=1).crossing(lambda x:x-100)
    with pytest.raises(RuntimeError):RootPolicy(maxiter=1).crossing(lambda x:x*x-2,initial=2)
    with pytest.raises(ValueError):RootPolicy().crossing(lambda x:np.nan)
    with pytest.raises(ValueError):RootPolicy().crossing(lambda x:x+1)
    with pytest.raises(ValueError):RootPolicy(rtol=1e-20)
    s=np.array([1.,2.]);b=np.array([4.,5.]);n=np.array([4,5])
    before=[x.copy() for x in (s,b,n)]
    result=PoissonCLsCounting().observed_upper_limit(s,b,n)
    for x,y in zip((s,b,n),before):np.testing.assert_array_equal(x,y)
    n[:]=100
    assert result.reference_counts==(4.,5.)
