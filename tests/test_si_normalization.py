"""Independent proton/coherent SI oracles; synthetic nuclear inputs only."""
from dataclasses import FrozenInstanceError, replace
from math import log, pi, sqrt

import numpy as np
import pytest
from scipy.integrate import quad

from wimp_legend.backgrounds import ZeroBackground
from wimp_legend.detector import ConstantQuenching, DetectorResponse, NoSmearing
from wimp_legend.interactions import CANONICAL_OPERATORS, NREFT_CONVENTION_ID
from wimp_legend.rates import NREFTRateConfig, differential_cross_section_GeV_minus3
from wimp_legend.sensitivity import SensitivityAnalysis
from wimp_legend.si import (SIO1Signal, si_proton_coefficient_GeV_minus2,
                            si_proton_cross_section_cm2, si_wilson_coefficients)
from wimp_legend.statistics import PoissonCLsCounting
from test_rates import synthetic, polynomial_pdf


def conversion():
    # SI defining constants, independent of the production conversion constant.
    hc = 6.62607015e-34*29979245800/(2*pi*1.602176634e-10)
    return hc*hc


def coherent_dataset(A=20, Z=4, mass=10.):
    # A,Z are synthetic charges, not imported physical isotope data.
    charges = (A, 2*Z-A)
    dataset = synthetic({('M', t, u): (2/(4*pi)*charges[t]*charges[u],)
                         for t in (0, 1) for u in (0, 1)})
    return replace(dataset, isotope=replace(dataset.isotope, A=A, Z=Z, mass_GeV=mass))


def signal(**changes):
    inputs = dict(sigma_p_ref_cm2=1e-40, m_chi_GeV=5., m_proton_GeV=1.,
        fn_over_fp=1., proton_sign=1, dataset=coherent_dataset(), rho_chi_GeV_cm3=.3,
        targets_per_kg=1e25, j_chi=.5, config=NREFTRateConfig(2.),
        speed_pdf=polynomial_pdf(.01), beta_max=.01)
    inputs.update(changes)
    return SIO1Signal(**inputs)


@pytest.mark.parametrize('m_chi,m_proton', [(5., 1.), (.2, 1.3), (1000., .8)])
@pytest.mark.parametrize('sign', [-1, 1])
def test_absolute_conversion_and_signed_roundtrip(m_chi, m_proton, sign):
    sigma = 2e-43
    mu = m_chi*m_proton/(m_chi+m_proton)
    expected = sign*sqrt(sigma*pi/(conversion()*mu*mu))
    actual = si_proton_coefficient_GeV_minus2(sigma, m_chi_GeV=m_chi,
        m_proton_GeV=m_proton, proton_sign=sign)
    assert actual == pytest.approx(expected, rel=8e-16)
    assert si_proton_cross_section_cm2(actual, m_chi_GeV=m_chi,
                                      m_proton_GeV=m_proton) == pytest.approx(sigma, rel=1e-15)


@pytest.mark.parametrize('ratio', [0., 1., -.5, 2.])
def test_proton_neutron_and_half_isospin_mapping(ratio):
    values = si_wilson_coefficients(1e-40, m_chi_GeV=5, m_proton_GeV=1,
                                    fn_over_fp=ratio, proton_sign=1)
    cp = sqrt(1e-40*pi/conversion())/(5/6)
    assert values.get(1, 0) == pytest.approx(cp*(1+ratio)/2, rel=1e-15, abs=0)
    assert values.get(1, 1) == pytest.approx(cp*(1-ratio)/2, rel=1e-15, abs=0)
    assert values.to_proton_neutron()[1] == pytest.approx((cp, ratio*cp), rel=1e-15, abs=0)
    assert all(values.get(op, tau) == 0 for op in CANONICAL_OPERATORS if op != 1 for tau in (0, 1))
    assert values.convention_id == NREFT_CONVENTION_ID


@pytest.mark.parametrize('scale', [1e-200, 1e200])
def test_conversion_avoids_intermediate_coefficient_square_overflow(scale):
    # Equal masses M and cp=1/M give mu*cp=1/2 for any representable M.
    actual = si_proton_cross_section_cm2(1/scale, m_chi_GeV=scale, m_proton_GeV=scale)
    assert actual == pytest.approx(conversion()/(4*pi), rel=8e-16)
    inverse = si_proton_coefficient_GeV_minus2(actual, m_chi_GeV=scale,
                                             m_proton_GeV=scale, proton_sign=1)
    assert inverse == pytest.approx(1/scale, rel=1e-15)


@pytest.mark.parametrize('ratio', [0., 1., -.5])
def test_free_proton_and_coherent_total_cross_section_absolute_normalization(ratio):
    model = signal(fn_over_fp=ratio)
    for dataset, A, Z in [(coherent_dataset(1, 1, mass=1.), 1, 1), (coherent_dataset(), 20, 4)]:
        target_mass = dataset.isotope.mass_GeV
        mu_target, mu_p = 5*target_mass/(5+target_mass), 5/6
        endpoint = 2*mu_target**2*.01**2/target_mass*1e6
        total = quad(lambda E: differential_cross_section_GeV_minus3(E, .01, m_chi_GeV=5,
            dataset=dataset, coefficients=model.coefficients, j_chi=.5, config=model.config)*1e-6,
            0, endpoint, epsabs=0, epsrel=1e-10)[0]*conversion()
        expected = 1e-40*(mu_target/mu_p)**2*(Z+ratio*(A-Z))**2
        assert total == pytest.approx(expected, rel=2e-14, abs=0)


@pytest.mark.parametrize('ratio', [0., 1., -.5])
def test_callable_signal_absolute_rate_shapes_and_cross_section_scaling(ratio):
    model = signal(fn_over_fp=ratio)
    energies = np.array([[0., 1., 10.], [40., 100., 300.]])
    mu_target, mu_p, length = 5*10/15, 5/6, .01
    endpoint = 2*mu_target**2*length**2/10*1e6
    eta = 3/(2*length)*np.maximum(1-energies/endpoint, 0)
    # Direct SI formula, substituting cp² = pi*sigma_p/(mu_p²*(hbar*c)²).
    expected = (1e25*.3/5*29979245800*86400*1e-6*10/(2*mu_p**2)
                *1e-40*(4+ratio*16)**2*eta)
    actual = model(energies)
    np.testing.assert_allclose(actual, expected, rtol=3e-12, atol=0)
    assert isinstance(model(1.), float)
    assert model(np.empty((0, 2))).shape == (0, 2)
    np.testing.assert_allclose(replace(model, sigma_p_ref_cm2=9e-40)(energies), 9*actual,
                               rtol=3e-12, atol=0)
    np.testing.assert_allclose(replace(model, proton_sign=-1)(energies), actual, rtol=3e-12, atol=0)
    np.testing.assert_allclose(replace(model, config=NREFTRateConfig(.7))(energies), actual,
                               rtol=3e-12, atol=0)


def test_sigma_limit_from_complete_generic_analysis_and_reference_invariance():
    analysis = SensitivityAnalysis(DetectorResponse(ConstantQuenching(1), NoSmearing()),
                                   ZeroBackground(), PoissonCLsCounting())
    cross_sections = []
    for reference in (1e-40, 4e-40):
        model = signal(sigma_p_ref_cm2=reference)
        result = analysis.observed_limit(model, [0, 1, 3], 2., (0, 3), [0, 0])
        mu_p, mu_target = 5/6, 5*10/15
        endpoint = 2*mu_target**2*.01**2/10*1e6
        prefactor = 1e25*.3/5*29979245800*86400*1e-6*10/(2*mu_p**2)*reference*20**2*3/(2*.01)
        primitive = lambda E: prefactor*(E-E*E/(2*endpoint))
        expected_counts = 2*np.diff([primitive(E) for E in [0, 1, 3]])
        np.testing.assert_allclose(result.signal_counts, expected_counts, rtol=3e-12)
        assert result.limit.mu_upper == pytest.approx(-log(.1)/sum(expected_counts), rel=3e-12)
        cross_sections.append(model.cross_section_limit_cm2(result.limit.mu_upper))
        assert cross_sections[-1] == pytest.approx(reference*result.limit.mu_upper, rel=1e-15)
        assert result.limit.method == 'poisson_cls_counting'
    assert cross_sections[0] == pytest.approx(cross_sections[1], rel=3e-12, abs=0)


def test_zero_couplings_and_no_sensitivity_limit():
    assert si_proton_coefficient_GeV_minus2(0, m_chi_GeV=5, m_proton_GeV=1, proton_sign=-1) == 0
    assert si_proton_cross_section_cm2(0, m_chi_GeV=5, m_proton_GeV=1) == 0
    model = signal(rho_chi_GeV_cm3=0)
    result = SensitivityAnalysis(DetectorResponse(ConstantQuenching(1), NoSmearing()),
        ZeroBackground(), PoissonCLsCounting()).expected_limit(model, [0, 1], 1, (0, 1))
    assert result.limit.status == 'no_sensitivity'
    assert model.cross_section_limit_cm2(result.limit.mu_upper) == np.inf


def test_convention_required_entries_and_immutability_are_preserved():
    model = signal(breakpoints=np.array([.003]))
    assert model.breakpoints == (.003,)
    assert model.config.m_N_reference_GeV == 2
    assert model.m_proton_GeV == 1 and model.convention_id == NREFT_CONVENTION_ID
    with pytest.raises(FrozenInstanceError):
        model.sigma_p_ref_cm2 = 1
    with pytest.raises(TypeError):
        model.coefficients.values[1] = (1, 0)
    bad = replace(model.dataset, metadata=replace(model.dataset.metadata, convention_id='different'))
    with pytest.raises(ValueError, match='convention_id'):
        signal(dataset=bad)
    entries = dict(model.dataset.responses)
    entries.pop(next(iter(entries)))
    with pytest.raises(KeyError):
        signal(dataset=replace(model.dataset, responses=entries))


@pytest.mark.parametrize('bad', [-1, np.nan, np.inf, True, '1', [1.]])
def test_invalid_si_inputs(bad):
    with pytest.raises(ValueError):
        si_proton_coefficient_GeV_minus2(bad, m_chi_GeV=5, m_proton_GeV=1, proton_sign=1)
    with pytest.raises(ValueError):
        signal(sigma_p_ref_cm2=bad)


@pytest.mark.parametrize('changes', [dict(m_chi_GeV=0), dict(m_proton_GeV=0), dict(fn_over_fp=np.inf),
    dict(proton_sign=0), dict(proton_sign=True), dict(proton_sign=1.), dict(sigma_p_ref_cm2=0),
    dict(beta_max=1), dict(j_chi=.3), dict(breakpoints=[.01]), dict(speed_pdf=3)])
def test_invalid_reference_configuration(changes):
    with pytest.raises((ValueError, TypeError)):
        signal(**changes)


def test_nonrepresentable_outputs_and_invalid_limits_fail_explicitly():
    with pytest.raises(ValueError, match='floating-point range'):
        si_proton_cross_section_cm2(1e300, m_chi_GeV=1, m_proton_GeV=1)
    with pytest.raises(ValueError, match='floating-point range'):
        si_proton_cross_section_cm2(1e-300, m_chi_GeV=1, m_proton_GeV=1)
    model = signal()
    for bad in [-1, -np.inf, np.nan, True, '1', [1.]]:
        with pytest.raises(ValueError):
            model.cross_section_limit_cm2(bad)
