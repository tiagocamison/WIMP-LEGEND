"""Work-in-progress low-energy acceptance regression tests for the poster workbench.

Run from repository root: python -m pytest -q tests/test_ionization_acceptance.py
"""
from copy import deepcopy
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import ndtr

from analyses.legend_poster import load_config, validate_config
from analyses.legend_poster.analysis import fold_component
from analyses.legend_poster.models import responses
from wimp_legend.backgrounds import TruthBackground
from wimp_legend.energy import EnergyDomain


def configured(mode='none', *, floor=0., e50=.020, width=.005):
    c = deepcopy(load_config())
    d = c['detector']['presets'][c['detector']['default']]
    d['ionization_acceptance'] = dict(model=mode, floor_keVee=floor,
                                      e50_keVee=e50, width_keVee=width)
    return validate_config(c)


def test_disabled_default_preserves_existing_response():
    original = configured('none')
    er, nr = responses(original)
    assert er.truth_efficiency is None and nr.truth_efficiency is None
    removed = deepcopy(original)
    del removed['detector']['presets']['benchmark']['ionization_acceptance']
    er2, nr2 = responses(validate_config(removed))
    edges = np.array([0.03, 0.1, 0.2])
    np.testing.assert_array_equal(er.bin_probabilities(.035, edges),
                                  er2.bin_probabilities(.035, edges))
    np.testing.assert_array_equal(nr.bin_probabilities(.3, edges),
                                  nr2.bin_probabilities(.3, edges))


def test_hard_cutoff_applied_in_visible_energy_for_er_and_nr():
    c = configured('hard')
    er, nr = responses(c)
    assert er.truth_efficiency(.019) == 0
    assert er.truth_efficiency(.020) == 1
    q = nr.quenching
    nuclear_low = .05
    assert q.mean_energy_keV(nuclear_low) < .020
    assert nr.truth_efficiency(nuclear_low) == 0
    # The mapped NR breakpoint has the requested visible energy.
    breakpoint, = nr.truth_efficiency.integration_breakpoints_keV(0., 1.)
    assert q.mean_energy_keV(breakpoint) == pytest.approx(.020, abs=1e-11)
    assert nr.truth_efficiency(breakpoint*1.01) == 1
    edges = np.array([.1, .2])
    assert np.all(er.expected_line_counts(.019, 1., edges, 1.) == 0)
    assert er.expected_line_counts(.030, 1., edges, 1.).sum() > 0


def test_erf_half_efficiency_and_strict_floor():
    c = configured('erf', floor=.020, e50=.030, width=.005)
    er, nr = responses(c)
    assert er.truth_efficiency(.019) == 0
    assert er.truth_efficiency(.030) == pytest.approx(.5)
    assert 0 < er.truth_efficiency(.025) < .5
    assert er.truth_efficiency(.035) > .5
    assert er.truth_efficiency.integration_breakpoints_keV(0.,.1)
    assert nr.truth_efficiency.integration_breakpoints_keV(0.,1.)


@pytest.mark.parametrize('invalid', [
    {'model': 'arctan', 'floor_keVee': 0., 'e50_keVee': .02, 'width_keVee': .005},
    {'model': 'erf', 'floor_keVee': .03, 'e50_keVee': .02, 'width_keVee': .005},
    {'model': 'erf', 'floor_keVee': 0., 'e50_keVee': .02, 'width_keVee': 0.},
    {'model': 'hard', 'floor_keVee': -1., 'e50_keVee': .02, 'width_keVee': .005},
    {'model': True, 'floor_keVee': 0., 'e50_keVee': .02, 'width_keVee': .005},
])
def test_invalid_acceptance_settings_rejected(invalid):
    c = configured('none')
    c['detector']['presets']['benchmark']['ionization_acceptance'] = invalid
    with pytest.raises(ValueError): validate_config(c)


@pytest.mark.parametrize('backend', ['reference', 'vectorized'])
def test_continuum_hard_cutoff_matches_independent_quadrature(backend):
    c = configured('hard', e50=.02)
    er, _ = responses(c)
    def flat_rate(e):
        return np.ones_like(np.asarray(e, dtype=float))
    component = TruthBackground(EnergyDomain.ER, er, flat_rate, (0., .08))
    bins = np.array([.1, .2])
    n = deepcopy(c['numerical'])
    n['folding_backend'] = backend
    observed = fold_component(component,bins,n,signal=True,refined=True)[0]
    def integrand(e):
        sigma = float(er.resolution.sigma_at(e))
        return .9 * (ndtr((.2-e)/sigma) - ndtr((.1-e)/sigma))
    expected = quad(integrand, .02, .08, epsabs=1e-12, epsrel=1e-10)[0]
    assert observed == pytest.approx(expected, rel=1e-3, abs=1e-9)


def test_erf_acceptance_decreases_counts_relative_to_no_cutoff():
    c_none = configured('none')
    c_erf = configured('erf', floor=.020, e50=.030, width=.005)
    bins = np.array([.1, .2])
    def flat_rate(e): return np.ones_like(np.asarray(e, dtype=float))
    counts = []
    for c in (c_none,c_erf):
        er,_ = responses(c)
        component = TruthBackground(EnergyDomain.ER, er, flat_rate, (0., .08))
        n = deepcopy(c['numerical']); n['folding_backend']='vectorized'
        counts.append(float(fold_component(component,bins,n,signal=True)[0]))
    assert 0 < counts[1] < counts[0]


@pytest.mark.parametrize('values', [dict(mode='none'), dict(width_keVee=0.),
    dict(e50_keVee=True), dict(floor_keVee=.1), dict(quenching=object())])
def test_direct_acceptance_construction_validates_inputs(values):
    from analyses.legend_poster.ionization import IonizationAcceptance
    params=dict(mode='erf',e50_keVee=.02,width_keVee=.005)
    params.update(values)
    with pytest.raises((ValueError,TypeError)):IonizationAcceptance(**params)


def test_rejected_truth_events_never_reach_kernel():
    from wimp_legend.detector import KernelResponse
    def forbidden(*args):raise AssertionError('rejected event entered smearing')
    response=KernelResponse(forbidden,includes_efficiency=False,
        truth_efficiency=lambda e:0.)
    np.testing.assert_array_equal(response.expected_line_counts(.01,1.,[0.,1.],1.),[0.])
    np.testing.assert_array_equal(response.expected_counts(lambda e:1.,[0.,1.],1.,(0.,.1)),[0.])


@pytest.mark.parametrize('backend',['reference','vectorized'])
def test_narrow_erf_turn_on_matches_independent_quadrature(backend):
    c=configured('erf',floor=.015,e50=.020,width=1e-6)
    er,_=responses(c);bins=np.array([.01,.03,.1])
    component=TruthBackground(EnergyDomain.ER,er,lambda e:np.ones_like(e),(0.,.08))
    n=deepcopy(c['numerical']);n['folding_backend']=backend
    actual=fold_component(component,bins,n,signal=True,refined=True)
    def integrand(e,j):
        sigma=float(er.resolution.sigma_at(e))
        probability=ndtr((bins[j+1]-e)/sigma)-ndtr((bins[j]-e)/sigma)
        return .9*ndtr((e-.020)/1e-6)*probability
    points=[.015,.019992,.019996,.020,.020004,.020008]
    expected=[quad(lambda e:integrand(e,j),.015,.08,points=points[1:],epsabs=1e-12,epsrel=1e-10)[0] for j in range(2)]
    np.testing.assert_allclose(actual,expected,rtol=1e-3,atol=1e-9)
