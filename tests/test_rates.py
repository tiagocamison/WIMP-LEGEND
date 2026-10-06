"""Analytic synthetic rate checks; no production nuclear coefficients/data.

Expected values derive directly from scalar Eq. (40)/(50) and elementary
integrals of 3 beta²/L³. They do not call production moments or R/W helpers.
"""
import ast
from dataclasses import FrozenInstanceError, replace
from decimal import Decimal, localcontext
from math import exp, pi, sqrt
from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import quad, IntegrationWarning

from wimp_legend import constants
from wimp_legend.halo import StandardHaloModel
from wimp_legend.interactions import NREFT_CONVENTION_ID, WilsonCoefficients
from wimp_legend.kinematics import minimum_speed_c, maximum_recoil_energy_keV
from wimp_legend.rates import (
    NREFTRateConfig, transition_probability, differential_cross_section_GeV_minus3,
    integrate_speed_flux, differential_rate_per_kg_day_keV,
)
from wimp_legend.targets import (
    SourceReference, Isotope, MomentumVariable, ResponseKey, ResponseMetadata,
    PolynomialResponse, NuclearResponseDataset,
)

# Independent explicit channel vocabulary (not imported from the kernel).
CHANNELS = ('M', 'Sigma_double_prime', 'Sigma_prime', 'Phi_double_prime',
            'Phi_double_prime_M', 'Phi_tilde_prime', 'Delta', 'Delta_Sigma_prime')


def synthetic(entries=None, *, spin=.5, momentum=None):
    source = SourceReference('Synthetic analytic rate fixture', 'tests/test_rates.py', 'synthetic')
    isotope = Isotope('X', 20, 4, 10, spin, dict.fromkeys(('identity', 'mass_GeV', 'spin'), source))
    metadata = ResponseMetadata('Anand Eq. (40), synthetic values', ('0', '1'),
        't0=1,t1=tau3; half-sum Wilson convention', 'Raw W; no outer Q',
        'dimensionless', momentum or MomentumVariable('q2_GeV2'), source, NREFT_CONVENTION_ID)
    responses = {ResponseKey(ch, t, u): PolynomialResponse((0,), 0, None, source)
                 for ch in CHANNELS for t in (0, 1) for u in (0, 1)}
    for key, coeffs in (entries or {}).items():
        responses[ResponseKey(*key)] = PolynomialResponse(tuple(coeffs), 0, None, source)
    return NuclearResponseDataset('synthetic-analytic-only', isotope, metadata, responses)


def polynomial_pdf(length):
    def pdf(beta):
        x = np.asarray(beta)/length
        return np.where((x >= 0) & (x <= 1), 3*x*x/length, 0.)
    return pdf


def physical_conversion_reference():
    # Direct SI definition calculation in high precision; deliberately do not
    # reuse any production constant/helper for the expected conversion.
    with localcontext() as ctx:
        ctx.prec = 60
        p = Decimal('3.14159265358979323846264338327950288419716939937510582097494')
        hc = Decimal('6.62607015e-34')*Decimal('29979245800')/(2*p*Decimal('1.602176634e-10'))
        return float(hc*hc*Decimal('1e-6')*Decimal('29979245800')*Decimal(86400))


def rate(energy, *, ds=None, raw=None, length=.01, m_chi=5, density=.3, count=1e25, config=None, **kw):
    return differential_rate_per_kg_day_keV(energy, dataset=ds if ds is not None else synthetic({('M', 0, 0): (1,)}),
        coefficients=WilsonCoefficients({1: (2e-7, 0)} if raw is None else raw),
        m_chi_GeV=m_chi, rho_chi_GeV_cm3=density, targets_per_kg=count,
        config=config or NREFTRateConfig(1), j_chi=.5,
        speed_pdf=polynomial_pdf(length), beta_max=length, **kw)


def test_si_conversion_from_defining_constants():
    assert constants.C_CM_S == 29979245800
    assert constants.SECONDS_PER_DAY == 86400
    assert constants.HBAR_C_GEV_CM == pytest.approx(1.973269804593025e-14, rel=5e-16)
    got = constants.GEV_MINUS2_TO_CM2*constants.KEV_TO_GEV*constants.C_CM_S*constants.SECONDS_PER_DAY
    assert got == pytest.approx(physical_conversion_reference(), rel=5e-16)
    assert constants.GEV_MINUS2_TO_CM2 == pytest.approx(3.8937937217185937e-28, rel=5e-16)


def test_config_required_immutable_and_explicit():
    with pytest.raises(TypeError):
        NREFTRateConfig()
    value = NREFTRateConfig(np.float64(1.5))
    assert type(value.m_N_reference_GeV) is float
    with pytest.raises(FrozenInstanceError):
        value.m_N_reference_GeV = 2
    # Config has no quadrature policy or separately invented convention ID.
    assert set(value.__dataclass_fields__) == {'m_N_reference_GeV'}


@pytest.mark.parametrize('bad', [0, -1, np.nan, np.inf, True, '1', [1]])
def test_config_invalid_mass(bad):
    with pytest.raises(ValueError):
        NREFTRateConfig(bad)


@pytest.mark.parametrize('spin', [0, .5, 1.5])
def test_o1_probability_and_cross_section(spin):
    # Sum = 4*2 + 6*(-3) + 6*5 + 9*7 = 83.
    ds = synthetic({('M', 0, 0): (2,), ('M', 0, 1): (-3,),
                    ('M', 1, 0): (5,), ('M', 1, 1): (7,)}, spin=spin)
    c = WilsonCoefficients({1: (2, 3)})
    expected = 4*pi/(2*spin+1)*83
    q2 = np.array([[0], [.02]])
    V = np.array([0., .001, .1])
    actual = transition_probability(ds, c, q2_GeV2=q2, V=V, j_chi=1, config=NREFTRateConfig(1))
    np.testing.assert_allclose(actual, np.full((2, 3), expected), rtol=5e-16)
    cross = differential_cross_section_GeV_minus3(2, .01, m_chi_GeV=5, dataset=ds,
             coefficients=c, j_chi=.5, config=NREFTRateConfig(1))
    assert cross == pytest.approx(10/(2*pi*.01**2)*expected, rel=5e-16)


@pytest.mark.parametrize('channel,raw,expected', [
    ('M', {1: (2, 0)}, lambda q,v: 4+0*q+0*v),
    ('Sigma_double_prime', {10: (2, 0)}, lambda q,v: q+0*v),
    ('Sigma_prime', {7: (2, 0)}, lambda q,v: v/2+0*q),
    ('Phi_double_prime', {3: (2, 0)}, lambda q,v: q*q+0*v),
    ('Phi_double_prime_M', {3: (2, 0), 1: (3, 0)}, lambda q,v: 6*q+0*v),
    ('Phi_tilde_prime', {12: (2, 0)}, lambda q,v: q/4+0*v),
    ('Delta', {8: (2, 0)}, lambda q,v: q+0*v),
    ('Delta_Sigma_prime', {5: (2, 0), 4: (3, 0)}, lambda q,v: 1.5*q+0*v),
])
def test_all_channels_outer_q_exactly_once(channel, raw, expected):
    ds = synthetic({(channel, 0, 0): (1,)})
    Q = np.array([[0.], [.03], [.7], [2.]])
    V = np.array([0., .002, .13])
    # Synthetic reference mass 2, so q²=4Q; J=1/2 -> nuclear factor 2pi.
    actual = transition_probability(ds, WilsonCoefficients(raw), q2_GeV2=4*Q, V=V,
                                    j_chi=.5, config=NREFTRateConfig(2))
    np.testing.assert_allclose(actual, 2*pi*expected(Q,V), rtol=6e-16, atol=0)


@pytest.mark.parametrize('channel,left,right,sign_factor', [
    ('Phi_double_prime_M', 3, 1, 1),
    ('Delta_Sigma_prime', 5, 4, .25),
    ('Delta_Sigma_prime', 8, 9, -.25),
])
def test_ordered_interference_contraction(channel, left, right, sign_factor):
    ds = synthetic({(channel, 0, 1): (7,), (channel, 1, 0): (-2,)})
    c = WilsonCoefficients({left: (2, 7), right: (3, 11)})
    actual = transition_probability(ds, c, q2_GeV2=.12, V=.03, j_chi=.5, config=NREFTRateConfig(2))
    assert actual == pytest.approx(2*pi*.03*sign_factor*(7*2*11-2*7*3), rel=5e-16)


@pytest.mark.parametrize('entry', [('M', 0, 0), ('Delta_Sigma_prime', 1, 0)])
def test_missing_entries_fail_even_for_zero_couplings(entry):
    ds = synthetic()
    responses = dict(ds.responses)
    del responses[ResponseKey(*entry)]
    ds = replace(ds, responses=responses)
    with pytest.raises(KeyError, match='missing required'):
        transition_probability(ds, WilsonCoefficients(), q2_GeV2=0, V=0, j_chi=0, config=NREFTRateConfig(1))
    with pytest.raises(KeyError):
        rate(1e9, ds=ds, raw={})


@pytest.mark.parametrize('bad_id', ['different', NREFT_CONVENTION_ID+' ', NREFT_CONVENTION_ID.upper()])
def test_exact_convention_gate_precedes_evaluation(bad_id, monkeypatch):
    ds = synthetic()
    ds = replace(ds, metadata=replace(ds.metadata, convention_id=bad_id))
    def forbidden(*args, **kwargs):
        raise AssertionError('nuclear evaluation must not happen')
    monkeypatch.setattr(NuclearResponseDataset, 'evaluate', forbidden)
    with pytest.raises(ValueError, match='convention_id'):
        transition_probability(ds, WilsonCoefficients(), q2_GeV2=-1, V=0,
                               j_chi=.5, config=NREFTRateConfig(1))
    with pytest.raises(ValueError, match='convention_id'):
        rate(np.empty((0, 2)), ds=ds)


@pytest.mark.parametrize('changes', [{'isospin_labels': ('p', 'n')}, {'response_units': 'fm^2'}])
def test_contradictory_metadata_rejected(changes):
    ds = synthetic()
    with pytest.raises(ValueError):
        rate(0, ds=replace(ds, metadata=replace(ds.metadata, **changes)))


def test_source_coordinate_and_validated_range_are_used():
    s = synthetic().metadata.source
    ds = synthetic({('M',0,0): (1,2)}, momentum=MomentumVariable('y', 4, s))
    c = WilsonCoefficients({1: (3, 0)})
    actual = transition_probability(ds, c, q2_GeV2=.25, V=0, j_chi=0, config=NREFTRateConfig(1))
    assert actual == pytest.approx(2*pi*9*3)  # y=1, W=1+2y
    responses = dict(ds.responses)
    key = ResponseKey('M', 0, 0)
    responses[key] = replace(responses[key], validated_range=(0, .5), validated_range_source=s)
    limited = replace(ds, responses=responses)
    with pytest.raises(ValueError, match='validated_range'):
        transition_probability(limited, c, q2_GeV2=.25, V=0, j_chi=0, config=NREFTRateConfig(1))


@pytest.mark.parametrize('reference_mass', [.5, 1, 2])
def test_q_dependent_o11_and_explicit_mass(reference_mass):
    ds = synthetic({('M',0,0): (1,)})
    c = WilsonCoefficients({11: (2, 0)})
    actual = transition_probability(ds, c, q2_GeV2=.04, V=0, j_chi=.5,
                                    config=NREFTRateConfig(reference_mass))
    assert actual == pytest.approx(2*pi*.04/reference_mass**2)
    assert ds.metadata.convention_id == NREFT_CONVENTION_ID


@pytest.mark.parametrize('operator', [1, 11, 8])
def test_closed_form_polynomial_pdf_rates(operator):
    L, m, target, c, rho, count = .01, 5., 10., 2e-7, .3, 1e25
    energy = np.array([[0., .5, 3.], [10., 40., 300.]])
    # Independent CM kinematics, no production helper in expected values.
    mu = m*target/(m+target)
    q2 = 2*target*energy*1e-6
    a2 = q2/(4*mu**2)
    t = np.minimum(a2/L**2, 1)
    eta = 3/(2*L)*(1-t)
    if operator == 1:
        integral = eta*c*c
    elif operator == 11:
        integral = eta*c*c*q2/4  # S/3=1/4 and reference mass 1
    else:
        # Integral f/beta * (beta²-a²) = 3 L/4 (1-a²/L²)^2.
        integral = c*c/4 * 3*L/4*(1-t)**2
    # W00=1, J=1/2 -> Ptot=2pi R, cancels 2pi in dsigma.
    expected = count*rho/m*physical_conversion_reference()*target*integral
    actual = rate(energy, raw={operator: (c, 0)})
    np.testing.assert_allclose(actual, expected, rtol=2e-12, atol=0)
    assert actual.shape == energy.shape


def test_o1_total_cross_section_coherent_normalization():
    # Purely synthetic W constructed from coherent charges, not nuclear data.
    A, Z, J, cp, cn = 20, 4, .5, 2e-7, -3e-7
    charges = (A, 2*Z-A)
    entries = {('M',t,u): ((2*J+1)/(4*pi)*charges[t]*charges[u],)
               for t in (0,1) for u in (0,1)}
    ds = synthetic(entries)
    c = WilsonCoefficients.from_proton_neutron({1:(cp,cn)})
    amplitude = cp*Z+cn*(A-Z)
    assert transition_probability(ds, c, q2_GeV2=0, V=0, j_chi=.5,
                                  config=NREFTRateConfig(1)) == pytest.approx(amplitude**2, rel=5e-16)
    mass, beta = 5., .01
    mu = mass*10/(mass+10)
    end_keV = 2*mu**2*beta**2/10*1e6
    total = quad(lambda E: differential_cross_section_GeV_minus3(E, beta, m_chi_GeV=mass,
                 dataset=ds, coefficients=c, j_chi=.5, config=NREFTRateConfig(1))*1e-6,
                 0, end_keV, epsabs=0, epsrel=1e-10)[0]
    assert total == pytest.approx(mu**2/pi*amplitude**2, rel=1e-14)


def test_kinematic_threshold_and_signed_results():
    ds = synthetic({('M',0,0): (1,)})
    args = dict(dataset=ds, coefficients=WilsonCoefficients({8:(2,0)}),
                config=NREFTRateConfig(1), m_chi_GeV=5, j_chi=.5)
    minimum = float(minimum_speed_c(3., 5., 10.))
    below, above = np.nextafter(minimum, 0), np.nextafter(minimum, 1)
    assert differential_cross_section_GeV_minus3(3., below, **args) == 0
    assert differential_cross_section_GeV_minus3(3., minimum, **args) == 0
    value = differential_cross_section_GeV_minus3(3., above, **args)
    expected = 10*(above-minimum)*(above+minimum)/above**2
    assert value > 0 and value == pytest.approx(expected, rel=5e-16)
    with pytest.raises(ValueError, match='undefined'):
        differential_cross_section_GeV_minus3(0, 0, **args)
    assert differential_cross_section_GeV_minus3(1, 0, **args) == 0
    neg = synthetic({('M',0,0): (-1,)})
    assert rate(1, ds=neg) == -rate(1)


def test_rate_endpoint_and_shapes():
    end = float(maximum_recoil_energy_keV(5, 10, .01))
    # Choose endpoint beta explicitly from the same elastic boundary, avoiding
    # an artificial mismatch between independently rounded inverse operations.
    bound = float(minimum_speed_c(end, 5, 10))
    assert rate(end, length=bound) == 0
    assert rate(end*1.001, length=bound) == 0
    assert rate(end*(1-1e-7), length=bound) > 0
    assert rate(np.empty((0,2))).shape == (0,2)
    assert isinstance(rate(1), float)
    ds = synthetic({('M',0,0): (1,)})
    e, b = np.array([[0,1,2.]]), np.array([[.001],[.01]])
    ec, bc = e.copy(), b.copy()
    args = dict(m_chi_GeV=5, dataset=ds, coefficients=WilsonCoefficients({1:(2,0)}),
                j_chi=.5, config=NREFTRateConfig(1))
    out = differential_cross_section_GeV_minus3(e,b,**args)
    assert out.shape == (2,3)
    for i,j in np.ndindex(out.shape):
        assert out[i,j] == differential_cross_section_GeV_minus3(float(e[0,j]),float(b[i,0]),**args)
    assert differential_cross_section_GeV_minus3(np.empty((0,2)),.01,**args).shape == (0,2)
    out[:] = -1
    np.testing.assert_array_equal(e,ec)
    np.testing.assert_array_equal(b,bc)
    assert transition_probability(ds,args['coefficients'],q2_GeV2=np.empty((0,1)),V=np.zeros(3),
                                  j_chi=.5,config=args['config']).shape == (0,3)


def test_density_and_target_count_are_explicit_linear_factors():
    assert rate(1, density=.6) == pytest.approx(2*rate(1))
    assert rate(1, count=3e25) == pytest.approx(3*rate(1))
    assert rate(1, density=0) == rate(1, count=0) == 0


def test_generic_integral_nonpolynomial_and_breakpoint():
    L = .007
    result = integrate_speed_flux(polynomial_pdf(L), lambda b: exp(b/L), beta_min=0, beta_max=L)
    assert result == pytest.approx(3*L*(6-2*exp(1)), rel=3e-15)
    piecewise = lambda b: 2/L if b < L/2 else 0
    result = integrate_speed_flux(piecewise, lambda b: 7., beta_min=0, beta_max=L, breakpoints=[L/2])
    assert result == pytest.approx(7*L/4, rel=5e-16)
    assert integrate_speed_flux(piecewise, lambda b: 7., beta_min=L, beta_max=L) == 0


def test_shm_callable_integration():
    h = StandardHaloModel.from_km_s()
    ds = synthetic({('M',0,0): (1,)})
    c = WilsonCoefficients({1:(2e-7,0)})
    # Independent zero-threshold eta from its standard erf expression.
    from math import erf
    z, y = h.vesc/h.v0, h.v_lab/h.v0
    norm = erf(z)-2*z*exp(-z*z)/sqrt(pi)
    eta = (erf(y)-2*y*exp(-z*z)/sqrt(pi))/(norm*h.v_lab)
    expected = 1e25*.3/5*physical_conversion_reference()*10*(2e-7)**2*eta
    actual = differential_rate_per_kg_day_keV(0, m_chi_GeV=5, rho_chi_GeV_cm3=.3,
        targets_per_kg=1e25, dataset=ds, coefficients=c, j_chi=.5, config=NREFTRateConfig(1),
        speed_pdf=h.lab_speed_pdf, beta_max=h.max_lab_speed, breakpoints=h.integration_breakpoints)
    assert actual == pytest.approx(expected, rel=2e-10)


@pytest.mark.parametrize('field', ['q2_GeV2','V'])
@pytest.mark.parametrize('bad', [-1,np.nan,np.inf,True,1j,'1'])
def test_invalid_transition_input(field,bad):
    kw=dict(q2_GeV2=.01,V=.001,j_chi=.5,config=NREFTRateConfig(1))
    kw[field]=bad
    with pytest.raises(ValueError):
        transition_probability(synthetic(),WilsonCoefficients(),**kw)


@pytest.mark.parametrize('field', ['E_nr_keV','m_chi_GeV','rho_chi_GeV_cm3','targets_per_kg','beta_max','j_chi'])
@pytest.mark.parametrize('bad', [-1,np.nan,np.inf,True])
def test_invalid_rate_input(field,bad):
    kw=dict(E_nr_keV=1,m_chi_GeV=5,rho_chi_GeV_cm3=.3,targets_per_kg=1e25,
            dataset=synthetic(),coefficients=WilsonCoefficients(),j_chi=.5,config=NREFTRateConfig(1),
            speed_pdf=polynomial_pdf(.01),beta_max=.01)
    kw[field]=bad
    with pytest.raises(ValueError):
        differential_rate_per_kg_day_keV(**kw)


@pytest.mark.parametrize('field,bad', [('beta',-1),('beta',np.inf),('beta',1),('beta',True),
                                      ('j_chi',.3),('m_chi_GeV',0),('E_nr_keV',np.nan)])
def test_invalid_cross_section(field,bad):
    kw=dict(E_nr_keV=1,beta=.01,m_chi_GeV=5,dataset=synthetic(),coefficients=WilsonCoefficients(),
            j_chi=.5,config=NREFTRateConfig(1))
    kw[field]=bad
    with pytest.raises(ValueError):
        differential_cross_section_GeV_minus3(**kw)


@pytest.mark.parametrize('field,bad', [('beta_min',-1),('beta_max',0),('beta_max',1),
    ('epsrel',0),('epsrel',np.nan),('epsabs',-1),('breakpoints',[0]),('breakpoints',[.01]),
    ('breakpoints',[np.nan]),('breakpoints',[[.001]])])
def test_invalid_quadrature_policy(field,bad):
    kw=dict(beta_min=0,beta_max=.01)
    kw[field]=bad
    with pytest.raises(ValueError):
        integrate_speed_flux(polynomial_pdf(.01),lambda b:1.,**kw)


@pytest.mark.parametrize('bad', [-1,np.nan,np.inf,True,[1.]])
def test_invalid_pdf_sample(bad):
    with pytest.raises(ValueError):
        integrate_speed_flux(lambda b:bad,lambda b:1.,beta_min=0,beta_max=.01)


@pytest.mark.parametrize('bad', [np.nan,np.inf,True,[1.],'1',1j])
def test_invalid_cross_section_sample(bad):
    with pytest.raises(ValueError):
        integrate_speed_flux(lambda b:100.,lambda b:bad,beta_min=0,beta_max=.01)


def test_divergent_quadrature_warns_as_error():
    with pytest.raises(IntegrationWarning):
        integrate_speed_flux(lambda b:100.,lambda b:1/b**2,beta_min=0,beta_max=.01)


def test_architectural_boundaries_and_generic_quadrature():
    root = Path(__file__).resolve().parents[1]/'src/wimp_legend'
    for folder in ('halo','interactions','targets'):
        for path in (root/folder).glob('*.py'):
            for node in ast.walk(ast.parse(path.read_text())):
                if isinstance(node,ast.ImportFrom):
                    assert 'rates' not in (node.module or '').split('.')
                elif isinstance(node,ast.Import):
                    assert all('rates' not in a.name.split('.') for a in node.names)
    source = (root/'rates/recoil.py').read_text()
    tree = ast.parse(source)
    for node in ast.walk(tree):
        if isinstance(node,ast.ImportFrom):
            assert not set((node.module or '').split('.')) & {'detector','backgrounds','statistics','sensitivity'}
            assert all(a.name not in ('velocity_moment','mean_inverse_speed','build_velocity_integral') for a in node.names)
    assert (root/'rates/kinematics.py').read_text() == ''
