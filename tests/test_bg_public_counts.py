"""Default-policy regression gates for CODE_AUDIT_007_BG_BENCHMARK.

No plotting fold is used. References integrate piecewise-linear beta data
against Gaussian CDFs or integrate the ES polynomial analytically.
"""
from dataclasses import replace
from math import exp, expm1, fsum, log, pi
from pathlib import Path
import json
import warnings

import numpy as np
import pytest
from numpy.polynomial.legendre import leggauss
from scipy.integrate import IntegrationWarning
from scipy.special import ndtr

from wimp_legend import _analysis
from wimp_legend._analysis import IntegrationPolicy
from wimp_legend.backgrounds import BenchmarkConfig, TabulatedBetaShape, build_benchmark
from wimp_legend.backgrounds.benchmark import benchmark_inputs
from wimp_legend.constants import GEV_MINUS2_TO_CM2, SECONDS_PER_DAY
from wimp_legend.detector import (
    ConstantEfficiency, DetectorResponse, ElectronRecoilResponse,
    GaussianResolution, LindhardQuenching, NoSmearing,
)


@pytest.mark.parametrize('limit', [1, 2, 7, 300])
def test_dense_breakpoints_preserve_piecewise_integral(limit):
    x = np.linspace(0, 1, 453)**1.3
    y = 1 + np.sin(67*x)**2
    expected = fsum(np.diff(x)*(y[:-1]+y[1:])/2)
    got = IntegrationPolicy(limit=limit).integrate(
        lambda e: np.interp(e, x, y), 0, 1,
        breakpoints=tuple(x[::-1])+tuple(x[50:70])+(-1, 2),
    )
    assert got == pytest.approx(expected, rel=3e-15)


def test_dense_breakpoint_calls_preserve_tolerance_budget(monkeypatch):
    calls = []
    original = _analysis.quad

    def record(callback, lower, upper, **kwargs):
        calls.append((lower, upper, kwargs.copy()))
        return original(callback, lower, upper, **kwargs)

    monkeypatch.setattr(_analysis, 'quad', record)
    policy = IntegrationPolicy(epsabs=1e-12, limit=7)
    x = np.linspace(0, 2, 50)**2
    assert policy.integrate(lambda t: t*t, 0, 4, breakpoints=x) == pytest.approx(64/3)
    assert len(calls) > 1
    assert calls[0][0] == 0 and calls[-1][1] == 4
    for (_, hi, _), (lo, _, _) in zip(calls, calls[1:]):
        assert hi == lo
    assert fsum(k['epsabs'] for _, _, k in calls) == pytest.approx(policy.epsabs, rel=3e-15)
    for lo, hi, kwargs in calls:
        assert kwargs['limit'] == policy.limit
        assert kwargs['epsrel'] == policy.epsrel
        assert kwargs['epsabs'] == pytest.approx(policy.epsabs*(hi-lo)/4)
        assert len(kwargs['points']) < policy.limit//2


def test_dense_breakpoints_still_raise_integration_warning(monkeypatch):
    def fail(*args, **kwargs):
        warnings.warn('deliberate quadrature failure', IntegrationWarning)
    monkeypatch.setattr(_analysis, 'quad', fail)
    with pytest.raises(IntegrationWarning, match='deliberate'):
        IntegrationPolicy(limit=3).integrate(lambda x: 1., 0, 1,
                                             breakpoints=np.linspace(0, 1, 15))


def test_tabulated_beta_interpolation_cache_cannot_be_mutated():
    source = TabulatedBetaShape([0, 1, 2], [1, 2, 0], 'synthetic')
    for grid in (source._energy_grid, source._rate_grid):
        with pytest.raises(ValueError):
            grid.setflags(write=True)
        with pytest.raises(ValueError):
            grid[0] = 99
    np.testing.assert_array_equal(source([0, .5, 1, 2]), [.4, .6, .8, 0])
    assert source == replace(source)


@pytest.fixture(scope='module')
def configured_models():
    record = json.loads((Path(__file__).resolve().parents[1]/
                        'examples/bg_benchmark_config.json').read_text())
    config = BenchmarkConfig(**record['scenarios']['legend_poster'])
    p = record['plotting']
    def sigma(e):
        return np.sqrt(p['noise_sigma_keV']**2 +
                       p['fano_factor']*p['pair_energy_keV']*np.asarray(e))
    efficiency = p['constant_detection_efficiency']
    er = ElectronRecoilResponse(GaussianResolution(sigma),
                               observed_efficiency=ConstantEfficiency(efficiency))
    nr = DetectorResponse(LindhardQuenching(32, p['lindhard_k']), GaussianResolution(sigma),
                          observed_efficiency=ConstantEfficiency(efficiency))
    direct = ElectronRecoilResponse(NoSmearing(),
                                   observed_efficiency=ConstantEfficiency(efficiency))
    assert er.policy == nr.policy == direct.policy == IntegrationPolicy()
    return (build_benchmark(config, er_response=er, nr_response=nr),
            build_benchmark(config, er_response=direct, nr_response=nr),
            config, sigma, efficiency)


def ar_mean_intensity(config):
    nuclear, _, _, _ = benchmark_inputs()
    decay = log(2)/nuclear['decays']['Ar39']['half_life_day']
    average = (exp(-decay*config.cooldown_day)*
               -expm1(-decay*config.run_duration_day)/(decay*config.run_duration_day))
    return (config.ar39_activity_Bq_per_kg_Ar*SECONDS_PER_DAY*
            config.argon_mass_per_kg_Ge*average*config.lar_transfer_probability*
            config.lar_veto_survival_probability)


def test_default_ar39_full_support_no_smearing(configured_models):
    _, direct, config, _, efficiency = configured_models
    component = direct.components['Ar39']
    _, _, ar, _ = benchmark_inputs()
    assert component.spectrum_breakpoints_keV == ar.energies_keV
    expected = ar_mean_intensity(config)*efficiency*365.25
    got = component.expected_counts([-1, 566], 365.25)
    np.testing.assert_allclose(got, [expected], rtol=2e-13)


def test_default_ar39_gaussian_bins_against_table_cdf_reference(configured_models):
    model, _, config, sigma, efficiency = configured_models
    _, _, ar, _ = benchmark_inputs()
    # Integrate every source segment on full support, not a truncated window.
    e = np.asarray(ar.energies_keV)
    density = np.asarray(ar.unnormalized_dndE)/ar.source_integral
    x, w = leggauss(12)
    nodes = ((e[:-1, None]+e[1:, None])/2 + np.diff(e)[:, None]*x/2).ravel()
    weights = (np.diff(e)[:, None]*w/2).ravel()
    bins = np.array([.1, 1., 15.])
    probability = np.diff(ndtr((bins-nodes[:, None])/sigma(nodes)[:, None]), axis=1)
    expected = 365.25*ar_mean_intensity(config)*efficiency*(
        (weights*np.interp(nodes, e, density)) @ probability)
    got = model.components['Ar39'].expected_counts(bins, 365.25)
    np.testing.assert_allclose(got, expected, rtol=2e-10, atol=1e-16)


def es_integral_reference(spec, lower_keV, upper_keV):
    """Integrate each flux node's polynomial above each occupied shell edge."""
    sw = spec.sin2_theta_W
    result = []
    for flux in spec.fluxes:
        energy = np.asarray(flux.energies_MeV)*1e-3
        pee = np.asarray(spec.survival_probability(flux.energies_MeV))
        gl = pee*(.5+sw)+(1-pee)*(-.5+sw)
        gl2 = pee*(.5+sw)**2+(1-pee)*(-.5+sw)**2
        a = gl2+sw**2
        b = -2*sw**2/energy-gl*sw*spec.m_e_GeV/energy**2
        c = sw**2/energy**2
        upper = np.minimum(upper_keV*1e-6, 2*energy**2/(spec.m_e_GeV+2*energy))
        for binding, occupancy in spec.electron_shells:
            lower = max(lower_keV, binding)*1e-6
            # Stable factored polynomial integral; zero for closed channels.
            dt = np.maximum(0., upper-lower)
            value = dt*(a+b*(upper+lower)/2+c*(upper**2+upper*lower+lower**2)/3)
            result.extend(occupancy*np.asarray(flux.weights_cm2_s)*value)
    return (fsum(result)*2*spec.G_F_GeV_minus2**2*spec.m_e_GeV/pi*
            GEV_MINUS2_TO_CM2*spec.targets_per_kg*SECONDS_PER_DAY)


def test_default_solar_es_no_smearing_against_analytic_integrals(configured_models):
    _, model, _, _, efficiency = configured_models
    component = model.components['solar_ES']
    spec = component.spectrum
    points = set(component.spectrum_breakpoints_keV)
    for flux in spec.fluxes:
        energy = np.asarray(flux.energies_MeV)*1e-3
        assert set(2*energy*energy/(spec.m_e_GeV+2*energy)*1e6) <= points
    expected = 365.25*efficiency*es_integral_reference(spec, .1, 15)
    np.testing.assert_allclose(component.expected_counts([.1, 15], 365.25),
                               [expected], rtol=2e-10)


def test_default_solar_es_gaussian_against_independent_audit(configured_models):
    model, _, _, _, _ = configured_models
    # CODE_AUDIT_007: raw piecewise-linear source flux integration, independent
    # of the production flux quadrature. Allow its documented representation error.
    got = model.components['solar_ES'].expected_counts([.1, 15], 365.25)
    np.testing.assert_allclose(got, [.0384948281], rtol=1e-6)


def test_default_benchmark_count_api_against_independent_audit(configured_models):
    model, _, _, _, _ = configured_models
    got = model.expected_counts([.1, 15], 365.25)
    np.testing.assert_allclose(got, [83.08333164], rtol=1e-7)


def test_dense_breakpoints_enforce_combined_error_estimate(monkeypatch):
    # Individually plausible local errors can exceed the requested global
    # tolerance after cancellation between signed pieces. Never accept that.
    answers = iter([(1., 5e-9), (-.9, 4e-9), (.01, 5e-11)])
    monkeypatch.setattr(_analysis, 'quad', lambda *args, **kwargs: next(answers))
    with pytest.raises(IntegrationWarning, match='global error budget'):
        IntegrationPolicy(limit=1).integrate(lambda x: 1., 0, 3, breakpoints=[1, 2])


def test_default_solar_cevens_gaussian_against_independent_audit(configured_models):
    model, _, _, _, _ = configured_models
    component = model.components['solar_CEvNS']
    points = set(component.spectrum_breakpoints_keV)
    for flux in component.spectrum.fluxes:
        e = np.asarray(flux.energies_MeV)*1e-3
        assert set(2*e*e/component.spectrum.m_target_GeV*1e6) <= points
    got = component.expected_counts([.1, 15], 365.25)
    np.testing.assert_allclose(got, [.417681325], rtol=1e-6)
