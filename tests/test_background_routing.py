"""Synthetic NR/ER/observed routing; elementary integral and point-measure oracles."""
from dataclasses import FrozenInstanceError
from math import erf, sqrt

import numpy as np
import pytest

from wimp_legend.backgrounds import (
    BinnedBackground, CallableBackground, CompositeBackground, FlatBackground,
    MonoenergeticLine, ObservedLineBackground, TruthBackground,
)
from wimp_legend.detector import (
    CallableEfficiency, ConstantEfficiency, ConstantQuenching, DetectorResponse,
    ElectronRecoilResponse, GaussianResolution, IntegrationPolicy, KernelResponse,
    MatrixResponse, NoSmearing,
)
from wimp_legend.energy import EnergyDomain
from wimp_legend.sensitivity import SensitivityAnalysis
from wimp_legend.statistics import PoissonCLsCounting


def line(energy, rate, domain):
    return MonoenergeticLine(energy, rate, domain)


def test_three_representations_and_per_component_responses_in_composite():
    nr = DetectorResponse(ConstantQuenching(.5), NoSmearing(),
        truth_efficiency=ConstantEfficiency(.5), observed_efficiency=ConstantEfficiency(.4))
    er = ElectronRecoilResponse(NoSmearing(), truth_efficiency=ConstantEfficiency(.25),
                               observed_efficiency=ConstantEfficiency(.8))
    nr_background = TruthBackground(EnergyDomain.NR, nr, spectrum=lambda E: 1.,
        spectrum_energy_range_keV=(0, 4), lines=[line(1.5, 3, EnergyDomain.NR)])
    er_background = TruthBackground(EnergyDomain.ER, er, spectrum=lambda E: 2.,
        spectrum_energy_range_keV=(0, 3), lines=[line(1.5, 5, EnergyDomain.ER)])
    observed = [FlatBackground(1), BinnedBackground([0, 1, 2, 3], [1, 2, 3],
        'fixed_expected_counts', 10), ObservedLineBackground([line(2, 2, EnergyDomain.OBSERVED)])]
    np.testing.assert_allclose(nr_background.expected_counts([0, 1, 2, 3], 10), [10, 4, 0], rtol=2e-14)
    np.testing.assert_allclose(er_background.expected_counts([0, 1, 2, 3], 10), [4, 14, 4], rtol=2e-14)
    combined = CompositeBackground([nr_background, er_background, *observed])
    np.testing.assert_allclose(combined.expected_counts([0, 1, 2, 3], 10), [25, 30, 37], rtol=2e-14)
    # The signal response has zero acceptance; it must not act on backgrounds.
    analysis = SensitivityAnalysis(DetectorResponse(ConstantQuenching(.1), NoSmearing(),
        observed_efficiency=ConstantEfficiency(0)), combined, PoissonCLsCounting())
    result = analysis.observed_limit(lambda E: 1., [0, 1, 2, 3], 10, (0, 1), [0, 0, 0])
    assert result.background_counts == pytest.approx([25, 30, 37])
    assert result.limit.status == 'no_sensitivity'


@pytest.mark.parametrize('domain', [EnergyDomain.NR, EnergyDomain.ER])
def test_gaussian_line_probabilities_and_efficiencies_without_truth_quadrature(domain, monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError('point lines must not use truth quadrature')
    monkeypatch.setattr(IntegrationPolicy, 'integrate', forbidden)
    options = dict(truth_efficiency=ConstantEfficiency(.5), observed_efficiency=ConstantEfficiency(.4))
    response = (DetectorResponse(ConstantQuenching(.5), GaussianResolution(.2), **options)
                if domain is EnergyDomain.NR else ElectronRecoilResponse(GaussianResolution(.2), **options))
    energy = 2. if domain is EnergyDomain.NR else 1.
    background = TruthBackground(domain, response, lines=[line(energy, 2., domain)])
    # Reconstructed mean=1, sigma=.2; X*line_rate*eps_truth*eps_obs=2.
    expected = 2*erf(1/sqrt(2))/2
    actual = background.expected_counts([.8, 1., 1.2], 5.)
    np.testing.assert_allclose(actual, [expected, expected], rtol=2e-14)
    assert sum(actual) < 2  # lost tails remain lost


def test_truth_and_observed_efficiency_coordinates_for_lines():
    nr = DetectorResponse(ConstantQuenching(.2), NoSmearing(),
        truth_efficiency=CallableEfficiency(lambda E_nr: E_nr/4),
        observed_efficiency=CallableEfficiency(lambda E_obs: E_obs))
    background = TruthBackground(EnergyDomain.NR, nr, lines=[line(2, 3, EnergyDomain.NR)])
    assert background.expected_counts([0, 1], 10)[0] == pytest.approx(6.)
    calls = []
    er = ElectronRecoilResponse(GaussianResolution(lambda E_er: calls.append(E_er) or .2))
    er.expected_line_counts(2, 1, [1.8, 2.2], 1)
    assert calls == [2.]
    assert not hasattr(er, 'quenching')
    with pytest.raises(TypeError):
        ElectronRecoilResponse(NoSmearing(), quenching=ConstantQuenching(.2))


def test_continuum_and_lines_have_separately_declared_support():
    response = ElectronRecoilResponse(NoSmearing())
    background = TruthBackground(EnergyDomain.ER, response, spectrum=lambda E: 2.,
        spectrum_energy_range_keV=(0, 1), lines=[line(2, 3, EnergyDomain.ER)])
    np.testing.assert_allclose(background.expected_counts([0, 1, 3], 4), [8, 12], rtol=2e-14)


def test_observed_lines_endpoints_and_direct_continuum_have_no_response_stage():
    lines = [line(E, rate, EnergyDomain.OBSERVED) for E, rate in
             [(-1, 100), (0, 1), (1, 2), (2, 3), (3, 100)]]
    np.testing.assert_array_equal(ObservedLineBackground(lines).expected_counts([0, 1, 2], 3), [3, 15])
    # A negative reconstructed point is allowed; a negative truth point is not.
    assert ObservedLineBackground([lines[0]]).expected_counts([-2, 0], 1)[0] == 100
    observed = CallableBackground(lambda E_obs: 1+E_obs, provenance='synthetic observed input')
    np.testing.assert_allclose(observed.expected_counts([0, 1, 2], 2), [3, 5], rtol=2e-14)


def test_domain_and_response_mismatches_are_rejected():
    nr = DetectorResponse(ConstantQuenching(.5), NoSmearing())
    er = ElectronRecoilResponse(NoSmearing())
    for domain, response in [(EnergyDomain.NR, er), (EnergyDomain.ER, nr),
                              (EnergyDomain.OBSERVED, er), ('keV_nr', nr), (True, nr)]:
        with pytest.raises(ValueError):
            TruthBackground(domain, response, spectrum=lambda E: 1., spectrum_energy_range_keV=(0, 1))
    with pytest.raises(ValueError):
        TruthBackground(EnergyDomain.NR, nr, lines=[line(1, 1, EnergyDomain.ER)])
    with pytest.raises(ValueError):
        ObservedLineBackground([line(1, 1, EnergyDomain.NR)])
    with pytest.raises(TypeError):
        MonoenergeticLine(1, 1, 'keV_nr')


def test_kernel_port_lines_and_matrix_point_rejection():
    response = KernelResponse(lambda E, bins: [.3, .4], includes_efficiency=True)
    background = TruthBackground(EnergyDomain.NR, response, lines=[line(1, 2, EnergyDomain.NR)])
    np.testing.assert_allclose(background.expected_counts([0, 1, 2], 5), [3, 4], rtol=2e-15)
    matrix = MatrixResponse([0, 1], [0, 1], [[.5]], False, 'synthetic uniform average')
    continuum = TruthBackground(EnergyDomain.NR, matrix, spectrum=lambda E: 2., spectrum_energy_range_keV=(0, 1))
    assert continuum.expected_counts([0, 1], 4)[0] == pytest.approx(4.)
    with pytest.raises(TypeError, match='point lines'):
        TruthBackground(EnergyDomain.NR, matrix, lines=[line(.5, 1, EnergyDomain.NR)])


def test_source_containers_are_immutable_and_outputs_owned():
    points = np.array([.4, .41])
    interval = np.array([0., 1.])
    lines = [line(.5, 2, EnergyDomain.ER)]
    response = ElectronRecoilResponse(NoSmearing(), observed_breakpoints=points)
    background = TruthBackground(EnergyDomain.ER, response, spectrum=lambda E: 1.,
        spectrum_energy_range_keV=interval, spectrum_breakpoints_keV=points, lines=lines)
    points[:] = 9; interval[:] = 9; lines.clear()
    assert response.observed_breakpoints == (.4, .41)
    assert background.spectrum_breakpoints_keV == (.4, .41)
    assert background.spectrum_energy_range_keV == (0, 1)
    result = background.expected_counts([0, 1], 1)
    assert result[0] == pytest.approx(3.) and result.flags.owndata
    result[:] = 9
    assert background.expected_counts([0, 1], 1)[0] == pytest.approx(3.)
    with pytest.raises(FrozenInstanceError):
        background.energy_domain = EnergyDomain.NR
    with pytest.raises(FrozenInstanceError):
        background.lines[0].energy_keV = 2


@pytest.mark.parametrize('bad', [-1, np.nan, np.inf, True, '1', [1.]])
def test_invalid_line_inputs(bad):
    with pytest.raises(ValueError):
        MonoenergeticLine(bad, 1, EnergyDomain.NR)
    with pytest.raises(ValueError):
        MonoenergeticLine(1, bad, EnergyDomain.ER)


def test_invalid_continuum_inputs_and_zero_exposure():
    response = ElectronRecoilResponse(NoSmearing())
    for kwargs in [{}, {'spectrum': 3, 'spectrum_energy_range_keV': (0, 1)},
                   {'spectrum': lambda E: 1}, {'spectrum': lambda E: 1, 'spectrum_energy_range_keV': (0, 1, 2)},
                   {'lines': [line(1, 1, EnergyDomain.ER)], 'spectrum_breakpoints_keV': [1]}]:
        with pytest.raises((ValueError, TypeError)):
            TruthBackground(EnergyDomain.ER, response, **kwargs)
    background = TruthBackground(EnergyDomain.ER, response, spectrum=lambda E: -1,
                                 spectrum_energy_range_keV=(0, 1))
    with pytest.raises(ValueError):
        background.expected_counts([0, 1], 1)
    np.testing.assert_array_equal(background.expected_counts([0, 1], 0), [0])
    with pytest.raises(ValueError):
        background.expected_counts([0, 1], -1)
