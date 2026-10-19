"""Physical nuclear inputs: source replay, normalization and convention checks.

All isotope masses here are synthetic test inputs, never production mass data.
Printed form factors are a separate, rounded primary-literature cross-check.
"""
from dataclasses import FrozenInstanceError, replace
import hashlib
import json
from math import exp, fsum, log10, floor, pi, sqrt
from pathlib import Path

import numpy as np
import pytest

from wimp_legend.targets import (
    Isotope, MomentumVariable, SourceReference, ResponseKey, TargetComposition,
    germanium_afh2014_momentum, load_germanium_afh2014,
)
from wimp_legend.targets import germanium as ge

ROOT = Path(ge.__file__).parent/'data/ge_afh2014_v1'
FIXTURES = Path(__file__).parent/'fixtures'
AS = (70, 72, 73, 74, 76)
CHANNELS = ('M', 'Sigma_double_prime', 'Sigma_prime', 'Phi_double_prime',
            'Phi_tilde_prime', 'Delta', 'Phi_double_prime_M', 'Delta_Sigma_prime')
SYNTHETIC = SourceReference('Synthetic test isotope properties', 'test_ge_nuclear_data.py', 'synthetic')


def isotope(A):
    return Isotope('Ge', A, 32, A*.9315, 4.5 if A == 73 else 0,
                   {k: SYNTHETIC for k in ('identity', 'mass_GeV', 'spin')})


def load(A):
    return load_germanium_afh2014(isotope(A), momentum=germanium_afh2014_momentum(A))


def raw(A):
    return json.loads((ROOT/f'ge{A}.json').read_text())


def matrix(dataset, channel, y):
    q2 = 4*y/dataset.metadata.momentum.b_GeV_inv**2
    return np.array([[dataset.evaluate(channel=channel, tau=t, tau_prime=u, q2_GeV2=q2)
                      for u in (0, 1)] for t in (0, 1)])


@pytest.mark.parametrize('A', AS)
def test_complete_ordered_data_hashes_domains_and_sources(A):
    dataset, data = load(A), raw(A)
    manifest_bytes = (ROOT/'manifest.json').read_bytes()
    manifest = json.loads(manifest_bytes)
    assert hashlib.sha256(manifest_bytes).hexdigest() == ge._MANIFEST_SHA256
    assert hashlib.sha256((ROOT/f'ge{A}.json').read_bytes()).hexdigest() == manifest['isotopes'][str(A)]['sha256']
    assert manifest['supplement_sha256'] == '943a403e216633840887bbfcbec55c56e0c5813b86335b49868acb89e74a49b9'
    assert len(dataset.responses) == 32
    assert list(dataset.responses) == [ResponseKey(ch, t, u) for ch in CHANNELS for t in (0, 1) for u in (0, 1)]
    assert dataset.metadata.convention_id == 'anand2014_prc89_065501_half_isospin_v1'
    assert dataset.metadata.response_units == 'dimensionless'
    assert data['primary_source']['density_row_count'] == (83 if A == 73 else 15)
    assert data['isotope']['mass_GeV'] is None
    for response in dataset.responses.values():
        assert response.domain == (0, None)
        assert response.validated_range is None
        assert response.exponential_decay == 2
        assert 'no particle-number rescaling' in response.source.notes
        assert 'source response' in response.source.locator
    assert json.loads(dataset.metadata.source.notes)['normalization']['particle_number_rescaling'] is False
    assert dataset.metadata.source.kind == 'implementation'
    assert 'pinned primary supplement' in dataset.metadata.source.upstream[0]


@pytest.mark.parametrize('A', AS)
def test_source_number_closure_and_coherent_o1(A):
    dataset, data = load(A), raw(A)
    T, J = A/2-32, isotope(A).spin
    proton_terms, neutron_terms = [], []
    # Independent q=0 M matrix element: sqrt(2j+1)/sqrt(4pi), diagonal K=0.
    for row in data['density_rows']:
        K, tau, Np, jp, N, j = row['quantum']
        if K != 0:
            continue
        assert (Np, jp) == (N, j)
        p = sqrt(2/(2*T+1))/2 if tau == 0 else -sqrt(6*T/((2*T+1)*(T+1)))/2
        term = float(row['psi_decimal'])*sqrt(j+1)/sqrt(2*J+1)
        proton_terms.append(term*p)
        neutron_terms.append(term*(p if tau == 0 else -p))
    Z_eff, N_eff = fsum(proton_terms), fsum(neutron_terms)
    assert Z_eff == pytest.approx(float(data['normalization']['effective_protons_q0']), rel=3e-15)
    assert N_eff == pytest.approx(float(data['normalization']['effective_neutrons_q0']), rel=3e-15)
    W = matrix(dataset, 'M', 0)
    for cp, cn in ((1,1), (1,0), (0,1), (1,-.7)):
        c = np.array([(cp+cn)/2, (cp-cn)/2])
        assert 4*pi/(2*J+1)*(c@W@c) == pytest.approx((cp*Z_eff+cn*N_eff)**2, rel=5e-15)
    # Detect accidental closure repair: the supplied source deficit stays visible.
    ratio = (Z_eff+N_eff)**2/A**2
    assert .999 < ratio < 1-1e-6
    assert 4*pi/(2*J+1)*W[0,0] != pytest.approx(A*A, rel=1e-6)


@pytest.mark.parametrize('A', AS)
def test_existing_rate_engine_coherent_o1_uses_unrescaled_source(A):
    from wimp_legend.interactions import WilsonCoefficients
    from wimp_legend.rates import NREFTRateConfig, transition_probability
    dataset, data = load(A), raw(A)
    cp, cn = 1e-7, -3e-8
    coefficients = WilsonCoefficients.from_proton_neutron({1:(cp,cn)})
    expected = (cp*float(data['normalization']['effective_protons_q0'])
                +cn*float(data['normalization']['effective_neutrons_q0']))**2
    actual = transition_probability(dataset,coefficients,q2_GeV2=0,V=.001,
                                    j_chi=.5,config=NREFTRateConfig(1.0))
    assert actual == pytest.approx(expected,rel=4e-15)


@pytest.mark.parametrize('A', (70, 72, 74, 76))
def test_even_isotope_exact_selection_rule_zeros(A):
    dataset = load(A)
    for ch in ('Sigma_prime', 'Sigma_double_prime', 'Phi_tilde_prime', 'Delta', 'Delta_Sigma_prime'):
        for t in (0, 1):
            for u in (0, 1):
                assert dataset.responses[ResponseKey(ch,t,u)].coefficients == (0.,)
                assert np.array_equal(dataset.evaluate(channel=ch,tau=t,tau_prime=u,q2_GeV2=[0,.01,.2]), [0,0,0])
    for ch in ('M', 'Phi_double_prime', 'Phi_double_prime_M'):
        assert any(dataset.responses[ResponseKey(ch,t,u)].coefficients != (0.,) for t in (0,1) for u in (0,1))


@pytest.mark.parametrize('A', AS)
def test_pn_to_half_isospin_factor_four_and_mixed_transpose(A):
    dataset, data = load(A), raw(A)
    S = np.array([[1,1],[1,-1]])
    y = .37
    for ch in CHANNELS:
        pn = np.zeros((2,2))
        for r in data['source_pn_polynomials']:
            if r['channel'] == ch:
                pn[r['nucleon'],r['nucleon_prime']] = np.polynomial.polynomial.polyval(y, list(map(float,r['coefficients'])))*exp(-2*y)
        source_W = (2*dataset.isotope.spin+1)/(16*pi)*S@pn@S.T
        if ch in ('Phi_double_prime_M', 'Delta_Sigma_prime'):
            source_W = source_W.T
        np.testing.assert_allclose(matrix(dataset,ch,y),4*source_W,rtol=3e-12,atol=1e-13)
    if A == 73:
        for ch in ('Phi_double_prime_M','Delta_Sigma_prime'):
            W = matrix(dataset,ch,y)
            assert not np.isclose(W[0,1],W[1,0],rtol=1e-3)


@pytest.mark.parametrize('A', AS)
def test_self_response_psd_vectorization_and_immutability(A):
    dataset = load(A)
    grid = np.array([[0,.1,.3],[.5,1,2]])
    q2 = 4*grid/dataset.metadata.momentum.b_GeV_inv**2
    before = q2.copy()
    for ch in CHANNELS[:6]:
        out = dataset.evaluate(channel=ch,tau=0,tau_prime=1,q2_GeV2=q2)
        assert out.shape == grid.shape
        for y in grid.flat:
            W = matrix(dataset,ch,float(y))
            np.testing.assert_allclose(W,W.T,rtol=0,atol=1e-13)
            assert np.linalg.eigvalsh(W).min() >= -1e-10
    np.testing.assert_array_equal(q2,before)
    with pytest.raises(FrozenInstanceError):
        dataset.metadata.convention_id = 'another'
    with pytest.raises(TypeError):
        dataset.responses[ResponseKey('M',0,0)] = None
    with pytest.raises(ValueError):
        dataset.evaluate(channel='M',tau=0,tau_prime=0,q2_GeV2=-1)
    with pytest.raises(ValueError):
        dataset.responses[ResponseKey('M',0,0)].evaluate(-.001)
    assert isinstance(dataset.evaluate(channel='M',tau=0,tau_prime=0,q2_GeV2=0),float)
    assert np.isfinite(matrix(dataset,'M',5)).all()  # no invented finite y_max


def test_mass_and_composition_independence_explicit_oscillator_choice():
    ge76 = isotope(76)
    momentum = germanium_afh2014_momentum(76)
    original = load_germanium_afh2014(ge76,momentum=momentum)
    another = load_germanium_afh2014(replace(ge76,mass_GeV=77),momentum=momentum)
    composition = TargetComposition(((ge76,.7),(isotope(70),.3)), 'number', SYNTHETIC)
    assert list(composition)[0][1] == .7
    assert original.isotope is ge76
    assert np.array_equal(matrix(original,'M',.25),matrix(another,'M',.25))
    b_fm = sqrt(41.467/(45*76**(-1/3)-25*76**(-2/3)))
    assert momentum.b_GeV_inv == b_fm*5.0677
    assert '1479' in momentum.b_source.locator
    alternate = MomentumVariable('y',10.,SYNTHETIC)
    changed = load_germanium_afh2014(ge76,momentum=alternate)
    assert changed.dataset_id != original.dataset_id
    assert changed.metadata.momentum is alternate
    np.testing.assert_allclose(matrix(changed,'M',.25),matrix(original,'M',.25),rtol=3e-15)
    with pytest.raises(TypeError):
        load_germanium_afh2014(ge76)  # mandatory oscillator input
    with pytest.raises(ValueError):
        load_germanium_afh2014(ge76,momentum=MomentumVariable('q2_GeV2'))
    with pytest.raises(ValueError):
        load_germanium_afh2014(replace(ge76,spin=1),momentum=momentum)
    with pytest.raises(ValueError):
        load_germanium_afh2014(replace(ge76,symbol='Zn',Z=30),momentum=momentum)
    for A in (71,True,76.0):
        with pytest.raises(ValueError):
            germanium_afh2014_momentum(A)


@pytest.mark.parametrize('filename', ('manifest.json', 'ge76.json'))
def test_loader_rejects_tampered_package_data(tmp_path, monkeypatch, filename):
    import shutil
    target = tmp_path/'data/ge_afh2014_v1'
    shutil.copytree(ROOT,target)
    with (target/filename).open('a') as stream:
        stream.write(' ')
    monkeypatch.setattr(ge,'files',lambda package: tmp_path)
    with pytest.raises(ValueError,match='checksum mismatch'):
        germanium_afh2014_momentum(76)


PRINTED = json.loads((FIXTURES/'ge_fitzpatrick_printed.json').read_text())
@pytest.mark.parametrize('record', PRINTED['records'], ids=lambda r: f"Ge{r['A']}_{r['channel']}_{r['nucleon']}{r['nucleon_prime']}")
def test_published_fitzpatrick_coefficients_with_rounding_intervals(record):
    data = raw(record['A'])
    extracted = next(r for r in data['source_pn_polynomials'] if (r['channel'],r['nucleon'],r['nucleon_prime']) == (record['channel'],record['nucleon'],record['nucleon_prime']))
    for k, token in enumerate(record['coefficients']):
        printed = float(token)
        # The appendix prints 2 significant digits, including trailing zeros
        # in large integers; omitted high-order powers remain in our export.
        half_last_digit = .5*10**(floor(log10(abs(printed)))-1)
        assert abs(float(extracted['coefficients'][k])-printed) <= half_last_digit+1e-11


@pytest.mark.parametrize('A', AS)
def test_original_supplement_interpreter_finite_momentum_reference(A):
    fixture = json.loads((FIXTURES/'ge_afh2014_source_replay.json').read_text())
    assert fixture['source_sha256'] == raw(A)['primary_source']['supplement_sha256']
    dataset = load(A)
    samples = [r for r in fixture['datasets'] if r['A']==A]
    assert {r['y'] for r in samples} == {.25,1.}
    for sample in samples:
        assert len(sample['responses']) == 32
        for r in sample['responses']:
            value = matrix(dataset,r['channel'],sample['y'])[r['tau'],r['tau_prime']]
            assert value == pytest.approx(r['value'],rel=4e-10,abs=2e-12)
