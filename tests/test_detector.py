"""Independent elementary integrals and probability references, synthetic only."""
from dataclasses import FrozenInstanceError, replace
from math import erf, sqrt, exp, pi
import numpy as np
import pytest
from scipy.integrate import quad, IntegrationWarning
from wimp_legend.detector import (
    IntegrationPolicy, ConstantQuenching, LindhardQuenching, NoSmearing,
    GaussianResolution, ConstantEfficiency, CallableEfficiency,
    DetectorResponse, KernelResponse, MatrixResponse,
)


def test_identity_nonuniform_bins_analytic_spectrum():
    d=DetectorResponse(ConstantQuenching(1),NoSmearing())
    bins=np.array([0.,.2,1.,3.])
    actual=d.expected_counts(lambda e:2+3*e,bins,7,(0,3))
    primitive=lambda e:2*e+1.5*e**2
    np.testing.assert_allclose(actual,7*np.diff(primitive(bins)),rtol=2e-14)
    assert sum(actual)==pytest.approx(7*(6+13.5))
    assert np.all(d.expected_counts(lambda e:1,bins,0,(0,3))==0)


@pytest.mark.parametrize('q',[1.,.2,.01])
def test_constant_quenching_jacobian_and_total(q):
    d=DetectorResponse(ConstantQuenching(q),NoSmearing())
    bins=np.array([0,.1,.3,1.])*q
    actual=d.expected_counts(lambda e:4,bins,2,(0,1))
    np.testing.assert_allclose(actual,8*np.diff(bins)/q,rtol=1e-14)
    assert sum(actual)==pytest.approx(8)


def test_lindhard_independent_formula_and_zero_shapes():
    e=np.array([[0.,.1],[1.,5.]])
    m=LindhardQuenching(17,.11)
    eps=11.5*e*17**(-7/3);g=3*eps**.15+.7*eps**.6+eps
    np.testing.assert_allclose(m.mean_energy_keV(e),e*.11*g/(1+.11*g),rtol=4e-16)
    assert m.mean_energy_keV(0)==0
    assert isinstance(m.mean_energy_keV(1),float)
    assert m.mean_energy_keV(np.empty((0,2))).shape==(0,2)
    assert np.all(np.diff(m.mean_energy_keV(np.linspace(0,10,100)))>0)
    d=DetectorResponse(m,NoSmearing())
    obs=m.mean_energy_keV([0,1,5])
    np.testing.assert_allclose(d.expected_counts(lambda e:2,obs,1,(0,5)),[2,8],rtol=2e-13)


@pytest.mark.parametrize('bad',[-1,np.nan,np.inf,True,1j,'1'])
def test_models_reject_invalid_energy(bad):
    for model in [ConstantQuenching(.2),LindhardQuenching(17,.11)]:
        with pytest.raises(ValueError):model.mean_energy_keV(bad)
    with pytest.raises(ValueError):GaussianResolution(1).bin_probabilities(bad,[0,1])


@pytest.mark.parametrize('bad',[0,-1,1.01,np.inf,np.nan,True])
def test_bad_quenching(bad):
    with pytest.raises(ValueError):ConstantQuenching(bad)


@pytest.mark.parametrize('args',[(0,.1),(1,0),(1,-1),(1.5,.1),(True,.1),(1,np.inf)])
def test_bad_lindhard(args):
    with pytest.raises(ValueError):LindhardQuenching(*args)


@pytest.mark.parametrize('model',[NoSmearing(),GaussianResolution(0),GaussianResolution(lambda mean:0)])
def test_delta_edge_assignment(model):
    np.testing.assert_array_equal(model.bin_probabilities([0,1,2,3],[0,1,2]),
                                  [[1,0],[0,1],[0,1],[0,0]])
    assert model.bin_probabilities(np.empty((0,2)),[0,1]).shape==(0,2,1)


def test_gaussian_probabilities_normalization_threshold_loss_and_tails():
    g=GaussianResolution(.3)
    bins=np.array([-3.,.4,1.,1.7,4.])
    actual=g.bin_probabilities(1,bins)
    cdf=lambda x:(1+erf((x-1)/(.3*sqrt(2))))/2
    np.testing.assert_allclose(actual,np.diff([cdf(x) for x in bins]),rtol=1e-12,atol=2e-16)
    assert sum(actual)==pytest.approx(1,abs=2e-14)
    assert g.bin_probabilities(1,[1,1.3])[0]==pytest.approx(erf(1/sqrt(2))/2)
    # Far positive tail independent PDF quadrature avoids subtracting near-one CDFs.
    tail=GaussianResolution(1).bin_probabilities(0,[9,9.2])[0]
    ref=quad(lambda z:exp(-z*z/2)/sqrt(2*pi),9,9.2,epsabs=0,epsrel=1e-12)[0]
    assert tail>0 and tail==pytest.approx(ref,rel=2e-14)
    assert sum(g.bin_probabilities(1,[.9,1.1]))<.3


@pytest.mark.parametrize('sigma',[1e-5,10.])
def test_narrow_and_broad_gaussian(sigma):
    actual=GaussianResolution(sigma).bin_probabilities(2,[2-sigma,2,2+sigma])
    np.testing.assert_allclose(actual,[erf(1/sqrt(2))/2]*2,rtol=2e-11)


def test_callable_sigma_means_and_ownership():
    g=GaussianResolution(lambda mean:.1+mean*.2)
    e=np.array([0.,1.,2.]);original=e.copy()
    d=DetectorResponse(ConstantQuenching(.2),g)
    result=d.bin_probabilities(e,[0,1,2])
    assert result.shape==(3,2)
    np.testing.assert_array_equal(e,original)
    assert g.sigma_at([0,.2,.4])==pytest.approx([.1,.14,.18])
    e[:]=9
    assert result[0,0]>.49
    with pytest.raises(FrozenInstanceError):g.sigma_keV=2


def test_truth_and_observed_efficiency_stages_and_in_bin_integration():
    base=DetectorResponse(ConstantQuenching(1),GaussianResolution(.2))
    eps=CallableEfficiency(lambda obs: .2 if obs<.4 else .8)
    observed=replace(base,observed_efficiency=eps,observed_breakpoints=[.4])
    mean=.5
    # Independent integral of piecewise efficiency against the normal PDF.
    density=lambda x:exp(-.5*((x-mean)/.2)**2)/(.2*sqrt(2*pi))
    expected=.2*quad(density,0,.4)[0]+.8*quad(density,.4,1)[0]
    assert observed.bin_probabilities(mean,[0,1])[0]==pytest.approx(expected,rel=2e-12)
    truth=replace(base,truth_efficiency=ConstantEfficiency(.4))
    np.testing.assert_allclose(truth.expected_counts(lambda e:1,[0,.5,1],2,(0,1)),
        .4*base.expected_counts(lambda e:1,[0,.5,1],2,(0,1)),rtol=1e-13)
    with pytest.raises(ValueError):CallableEfficiency(lambda e:1.01)(1)
    assert CallableEfficiency(lambda e:.2 if e<1 else .4)([0,1]).tolist()==[.2,.4]
    assert CallableEfficiency(lambda e:.5)([]).shape==(0,)


def test_kernel_matrix_equivalence_and_efficiency_ownership():
    matrix=np.array([[.8,.2],[.1,.6]])
    d=MatrixResponse([0,1,2],[0,2,3],matrix,False,'piecewise constant probabilities')
    kernel=KernelResponse(lambda e,edges: matrix[:,0 if e<1 else 1],False)
    truth=lambda e:1+e
    actual=d.expected_counts(truth,[0,2,3],2,(0,2))
    ref=kernel.expected_counts(truth,[0,2,3],2,(0,2),truth_breakpoints=[1])
    np.testing.assert_allclose(actual,ref,rtol=2e-14)
    np.testing.assert_allclose(actual,matrix@np.array([3.,5.]),rtol=1e-14)
    matrix[:]=0
    assert sum(d.apply_truth_counts([3,5]))>0
    with pytest.raises(FrozenInstanceError):d.matrix=()
    with pytest.raises(ValueError):replace(d,includes_efficiency=True,truth_efficiency=ConstantEfficiency(.5))
    with pytest.raises(ValueError):KernelResponse(lambda e,edges:[1,0],True,ConstantEfficiency(.5))
    with pytest.raises(ValueError):d.expected_counts(truth,[0,1,3],1,(0,2))
    with pytest.raises(ValueError):d.expected_counts(truth,[0,2,3],1,(0,1))


@pytest.mark.parametrize('matrix',[[[1.1]], [[-1]], [[np.nan]], [[.6],[.6]], [[1,0]]])
def test_invalid_matrix(matrix):
    with pytest.raises(ValueError):MatrixResponse([0,1],[0,1],matrix,False,'test')


@pytest.mark.parametrize('bad',[[0,0],[1,0],[0,np.inf],[0,np.nan],[0],[True,False]])
def test_invalid_edges(bad):
    with pytest.raises(ValueError):DetectorResponse(ConstantQuenching(1),NoSmearing()).expected_counts(
        lambda e:1,bad,1,(0,2))


@pytest.mark.parametrize('bad',[-1,np.nan,np.inf,True,[1,2]])
def test_invalid_spectrum_sample(bad):
    with pytest.raises(ValueError):DetectorResponse(ConstantQuenching(1),NoSmearing()).expected_counts(
        lambda e:bad,[0,1],1,(0,1))


@pytest.mark.parametrize('probability',[[1.1],[-1],[np.nan],[.6,.6],[]])
def test_invalid_custom_kernel(probability):
    with pytest.raises(ValueError):KernelResponse(lambda e,edges:probability,False).expected_counts(
        lambda e:1,[0,1],1,(0,1))


def test_policy_warning_propagation(monkeypatch):
    import warnings
    import wimp_legend._analysis as module
    def warning(*args,**kwargs):
        warnings.warn('injected integration failure',IntegrationWarning)
        return (1.,0.)
    monkeypatch.setattr(module,'quad',warning)
    with pytest.raises(IntegrationWarning):IntegrationPolicy().integrate(lambda e:1,0,1)
    with pytest.raises(ValueError):IntegrationPolicy(epsrel=0)


def test_energy_dependent_efficiencies_use_the_declared_stages():
    detector=DetectorResponse(ConstantQuenching(.5),NoSmearing(),
        truth_efficiency=CallableEfficiency(lambda E_nr:E_nr/2),
        observed_efficiency=CallableEfficiency(lambda E_obs:E_obs))
    # epsilon_truth*epsilon_obs=E_nr^2/4, edges map to E_nr=[0,.6,2].
    expected=np.diff(np.array([0.,.6,2.])**3/12)
    np.testing.assert_allclose(detector.expected_counts(lambda e:1,[0,.3,1],1,(0,2)),
                               expected,rtol=2e-14)


def test_weighted_gaussian_positive_tail_stays_nonzero():
    expected=.5*quad(lambda z:exp(-z*z/2)/sqrt(2*pi),9,9.2,epsabs=0,epsrel=1e-12)[0]
    actual=GaussianResolution(1).bin_probabilities(0,[9,9.2],
        observed_efficiency=CallableEfficiency(lambda obs:.5))[0]
    assert actual>0 and actual==pytest.approx(expected,rel=3e-14)


@pytest.mark.parametrize('bad',[-1,np.nan,np.inf,True,'1',1j])
def test_invalid_resolution_and_efficiency_parameters(bad):
    with pytest.raises(ValueError):GaussianResolution(bad)
    with pytest.raises(ValueError):ConstantEfficiency(bad)
    with pytest.raises(ValueError):GaussianResolution(lambda mean:bad).bin_probabilities(1,[0,2])
