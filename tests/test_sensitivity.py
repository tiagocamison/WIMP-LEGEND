"""Independent synthetic end-to-end analysis references; no real detector data."""
from dataclasses import FrozenInstanceError
from math import erf, exp, log, pi, sqrt
import ast
from pathlib import Path
import numpy as np
import pytest
from scipy.optimize import brentq
from scipy.stats import norm
from wimp_legend.detector import (
    DetectorResponse, ConstantQuenching, NoSmearing, GaussianResolution,
    ConstantEfficiency, CallableEfficiency, KernelResponse,
)
from wimp_legend.backgrounds import ZeroBackground, FlatBackground, CallableBackground
from wimp_legend.statistics import PoissonCLsCounting, BinnedPoissonProfile, LimitResult
from wimp_legend.sensitivity import SensitivityAnalysis, scan_parameter, ParameterScan


def normal_cdf(z):return .5*(1+erf(z/sqrt(2)))
def normal_primitive(z):return z*normal_cdf(z)+exp(-z*z/2)/sqrt(2*pi)


def flat_gaussian_counts(edges,*,q,sigma,length,rate,exposure,efficiency):
    # Integrate the normal CDF analytically over a flat truth interval [0,length].
    integral=lambda edge:sigma/q*(normal_primitive(edge/sigma)-normal_primitive((edge-q*length)/sigma))
    return np.diff([integral(edge) for edge in edges])*rate*exposure*efficiency


def cdf_by_recursion(n,mean):
    term=exp(-mean);result=term
    for k in range(1,n+1):term*=mean/k;result+=term
    return result


def test_identity_zero_background_observed_and_expected():
    analysis=SensitivityAnalysis(DetectorResponse(ConstantQuenching(1),NoSmearing()),
                                 ZeroBackground(),PoissonCLsCounting())
    expected=analysis.expected_limit(lambda e:3,[0,.2,1,2],5,(0,2))
    np.testing.assert_allclose(expected.signal_counts,[3,12,15],rtol=2e-14)
    assert expected.limit.mu_upper==pytest.approx(-log(.1)/30)
    observed=analysis.observed_limit(lambda e:3,[0,.2,1,2],5,(0,2),[0,0,0])
    assert observed.limit.mu_upper==expected.limit.mu_upper
    assert analysis.expected_limit(lambda e:3,[0,2],0,(0,2)).limit.status=='no_sensitivity'
    with pytest.raises(FrozenInstanceError):expected.signal_counts=()


def test_full_quenching_gaussian_efficiency_background_pipeline():
    bins=[0.,.2,.6,1.2];exposure=100.
    detector=DetectorResponse(ConstantQuenching(.2),GaussianResolution(.1),
        truth_efficiency=ConstantEfficiency(.4),observed_efficiency=ConstantEfficiency(.6))
    expected_signal=flat_gaussian_counts(bins,q=.2,sigma=.1,length=5,rate=2,
                                         exposure=exposure,efficiency=.24)
    expected_background=3*exposure*np.diff(bins)
    for backend in [PoissonCLsCounting(),BinnedPoissonProfile()]:
        result=SensitivityAnalysis(detector,FlatBackground(3),backend).expected_limit(
            lambda e:2,bins,exposure,(0,5))
        np.testing.assert_allclose(result.signal_counts,expected_signal,rtol=3e-12)
        np.testing.assert_allclose(result.background_counts,expected_background,rtol=3e-15)
        if isinstance(backend,PoissonCLsCounting):
            bg=sum(expected_background)
            n=0
            while cdf_by_recursion(n,bg)<.5:n+=1
            oracle=brentq(lambda m:cdf_by_recursion(n,bg+m*sum(expected_signal))/
                         cdf_by_recursion(n,bg)-.1,0,5)
            assert result.limit.kind=='expected_median'
        else:
            oracle=brentq(lambda m:2*sum(m*s-b*log(1+m*s/b) for s,b in
                        zip(expected_signal,expected_background))-norm.ppf(.9)**2,1e-5,5)
            assert result.limit.kind=='expected_asimov'
        assert result.limit.mu_upper==pytest.approx(oracle,rel=2e-10)


def test_narrow_truth_gaussian_migration_and_partial_window():
    mean,halfwidth,sigma=1.2,1e-6,.2
    d=DetectorResponse(ConstantQuenching(1),GaussianResolution(sigma))
    bins=[-1.,1.1,1.2,1.6,4.]
    interval=(mean-halfwidth,mean+halfwidth)
    # Normalize to the actual floating endpoints, not their nominal difference.
    density=1/(interval[1]-interval[0])
    actual=d.expected_counts(lambda e:density,bins,1,interval)
    expected=np.diff([normal_cdf((edge-mean)/sigma) for edge in bins])
    # Finite uniform source vs delta limit: error bounded at O((halfwidth/sigma)^2).
    np.testing.assert_allclose(actual,expected,atol=.2*(halfwidth/sigma)**2,rtol=0)
    partial=d.expected_counts(lambda e:density,[1.1,1.2],1,interval)
    assert sum(partial)<.2


def test_narrow_migration_window_discovered_by_known_gaussian_breakpoints():
    d=DetectorResponse(ConstantQuenching(1),GaussianResolution(1e-5))
    actual=d.expected_counts(lambda e:1,[1.9,1.90001],1,(0,10))
    assert actual[0]==pytest.approx(.00001,rel=3e-10)


def test_nonconstant_observed_efficiency_full_fold_independent_integral():
    from scipy.integrate import quad
    detector=DetectorResponse(ConstantQuenching(.5),GaussianResolution(.2),
        observed_efficiency=CallableEfficiency(lambda obs: .3 if obs<.4 else .7),
        observed_breakpoints=[.4])
    # Swap the two integrations. Truth E in [0,2], observed density is a difference
    # of normal CDFs divided by q. Integrate efficiency in observed energy.
    density=lambda obs:2*(normal_cdf(obs/.2)-normal_cdf((obs-1)/.2))
    expected=3*(.3*quad(density,0,.4,epsabs=1e-12)[0]+
                .7*quad(density,.4,1.3,epsabs=1e-12)[0])
    actual=detector.expected_counts(lambda e:1,[0,1.3],3,(0,2))
    assert actual[0]==pytest.approx(expected,rel=3e-11)


def test_external_backends_and_custom_kernel_grid():
    class CollaborationBackend:
        def observed_upper_limit(self,s,b,n,*,confidence_level):
            return LimitResult(7.,'synthetic_external','test-only construction',confidence_level,
                               'observed','test-only','success',True,tuple(n))
        def expected_upper_limit(self,s,b,*,confidence_level):
            return LimitResult(11.,'synthetic_external','test-only construction',confidence_level,
                               'expected','test-only','success',True)
        def acceptance_margin(self,s,b,n,*,confidence_level):return 1-sum(s)
        def acceptance_metadata(self):return ('synthetic_external','1-sum(S)','test-only')
    kernel=KernelResponse(lambda e,bins:[.2,.6],True)
    analysis=SensitivityAnalysis(kernel,CallableBackground(lambda e:1),CollaborationBackend())
    grid=analysis.expected_grid({'one':lambda e:1,'two':lambda e:2},[0,1,3],[2,4],(0,1))
    assert len(grid)==4 and all(point.analysis.limit.mu_upper==11 for point in grid)
    np.testing.assert_allclose(grid[0].analysis.signal_counts,[.4,1.2],rtol=1e-14)
    assert analysis.observed_limit(lambda e:1,[0,1,3],2,(0,1),[0,0]).limit.mu_upper==7
    scan=scan_parameter(lambda t:[t*t],parameter_domain=(-2,2),parameter_grid=[-2,-1,0,1,2],
        background_counts=[0],observed_counts=[0],backend=CollaborationBackend())
    assert scan.accepted==(False,True,True,True,False)
    assert scan.method=='synthetic_external'


def test_nonmonotonic_parameter_scan_no_unique_upper_limit():
    grid=np.linspace(-3,3,61)
    result=scan_parameter(lambda t:[5*(t*t-1)**2],parameter_domain=(-3,3),parameter_grid=grid,
        background_counts=[0],observed_counts=[0],backend=PoissonCLsCounting())
    expected=tuple(5*(t*t-1)**2<=-log(.1) for t in grid)
    assert result.accepted==expected and len(result.accepted_sample_ranges)==2
    assert result.status=='completed_explicit_grid' and 'not continuous' in result.notes[1]
    saved=result.parameters;grid[:]=0;assert result.parameters==saved
    with pytest.raises(FrozenInstanceError):result.parameter_domain=(0,1)


@pytest.mark.parametrize('kwargs',[
    {'parameter_domain':(0,2),'parameter_grid':[0,1]},
    {'parameter_grid':[0,0,1]}, {'parameter_domain':(0,0)},
    {'observed_counts':[.5]}, {'background_counts':[np.nan]},
    {'signal_counts':lambda t:[-1]},
])
def test_invalid_parameter_scan(kwargs):
    inputs=dict(signal_counts=lambda t:[t],parameter_domain=(0,1),parameter_grid=[0,.5,1],
                background_counts=[0],observed_counts=[0],backend=PoissonCLsCounting())
    inputs.update(kwargs)
    with pytest.raises(ValueError):scan_parameter(**inputs)


def test_existing_nreft_rate_as_supplied_truth_adapter():
    from test_rates import synthetic,polynomial_pdf
    from wimp_legend.rates import differential_rate_per_kg_day_keV,NREFTRateConfig
    from wimp_legend.interactions import WilsonCoefficients
    dataset=synthetic({('M',0,0):(1,)})
    rate=lambda e:differential_rate_per_kg_day_keV(e,m_chi_GeV=5,rho_chi_GeV_cm3=.3,
        targets_per_kg=1e25,dataset=dataset,coefficients=WilsonCoefficients({1:(2e-7,0)}),
        j_chi=.5,config=NREFTRateConfig(1),speed_pdf=polynomial_pdf(.01),beta_max=.01)
    detector=DetectorResponse(ConstantQuenching(1),NoSmearing())
    actual=detector.expected_counts(rate,[0,1,3],2,(0,3))
    # Independent SI + eta(Er)=3/(2L)*(1-Er/Emax); constant W and J=.5.
    mu=5*10/(5+10);endpoint=2*mu**2*.01**2/10*1e6
    hc=6.62607015e-34*29979245800/(2*pi*1.602176634e-10)
    prefactor=1e25*.3/5*hc**2*1e-6*29979245800*86400*10*(2e-7)**2*3/(2*.01)
    primitive=lambda e:prefactor*(e-e*e/(2*endpoint))
    np.testing.assert_allclose(actual,2*np.diff([primitive(e) for e in [0,1,3]]),rtol=4e-12)


def test_new_dependency_boundaries():
    root=Path(__file__).resolve().parents[1]/'src/wimp_legend'
    for area in ('detector','backgrounds','statistics'):
        forbidden={'rates','halo','targets','interactions','kinematics','sensitivity'}
        if area!='detector':forbidden.add('detector')
        if area!='backgrounds':forbidden.add('backgrounds')
        if area!='statistics':forbidden.add('statistics')
        for path in (root/area).glob('*.py'):
            for node in ast.walk(ast.parse(path.read_text())):
                names=[x.name for x in node.names] if isinstance(node,ast.Import) else [node.module or ''] if isinstance(node,ast.ImportFrom) else []
                for name in names:assert not set(name.split('.'))&forbidden,(path,name)
    for area in ('rates','interactions'):
        for path in (root/area).glob('*.py'):
            for node in ast.walk(ast.parse(path.read_text())):
                if isinstance(node,ast.ImportFrom):
                    assert not set((node.module or '').split('.'))&{'detector','backgrounds','statistics','sensitivity'}
