"""GA-AUDIT-001/002 regressions with independent synthetic-only references."""
from dataclasses import replace
from math import erfc, erf, exp, log, sqrt
from statistics import NormalDist

import numpy as np
import pytest
from scipy.optimize import brentq

from wimp_legend.backgrounds import FlatBackground
from wimp_legend.detector import (
    CallableEfficiency, ConstantQuenching, DetectorResponse, GaussianResolution,
    IntegrationPolicy, LindhardQuenching, NoSmearing,
)
from wimp_legend.sensitivity import SensitivityAnalysis
from wimp_legend.statistics import BinnedPoissonProfile, PoissonCLsCounting, RootPolicy


def normal_interval(lower, upper):
    """Independent erfc reference; never subtract CDFs close to one."""
    if lower >= 0:
        return .5*(erfc(lower/sqrt(2))-erfc(upper/sqrt(2)))
    if upper <= 0:
        return .5*(erfc(-upper/sqrt(2))-erfc(-lower/sqrt(2)))
    return .5*(erf(upper/sqrt(2))-erf(lower/sqrt(2)))


@pytest.mark.parametrize('sigma', [None, 0., 1e-5])
def test_GA001_declared_observed_window_not_lost(sigma):
    resolution = NoSmearing() if sigma is None else GaussianResolution(sigma)
    detector = DetectorResponse(ConstantQuenching(1), resolution,
        observed_efficiency=CallableEfficiency(lambda y: 1. if .4 < y < .41 else 0.),
        observed_breakpoints=(.4, .41))
    assert detector.expected_counts(lambda E: 1., [0, 1], 1., (0, 1))[0] == pytest.approx(
        .01, rel=1e-9, abs=0)


@pytest.mark.parametrize('mean', [.30, .32, .33])
@pytest.mark.parametrize('side', ['right', 'left'])
def test_GA002_weighted_tail_partition_invariance(mean, side):
    if side == 'left':
        mean = .81-mean
    samples = []

    def window(y):
        assert np.isfinite(y) and -10 <= y <= 10
        samples.append(y)
        return 1. if .4 < y < .41 else 0.

    efficiency = CallableEfficiency(window)
    reference = normal_interval((.4-mean)/.01, (.41-mean)/.01)
    assert reference > 0
    model = GaussianResolution(.01)
    fine = model.bin_probabilities(mean, [.4, .41], observed_efficiency=efficiency)[0]
    assert fine == pytest.approx(reference, rel=2e-12, abs=0)
    for bins in ([0, 1], [-10, 10]):
        coarse = model.bin_probabilities(mean, bins, observed_efficiency=efficiency,
                                         observed_breakpoints=[.4, .41])[0]
        assert coarse == pytest.approx(reference, rel=1e-7, abs=0)
    assert samples


def test_GA002_full_fold_with_explicit_truth_breakpoints():
    detector = DetectorResponse(ConstantQuenching(1), GaussianResolution(.01),
        observed_efficiency=CallableEfficiency(lambda y: 1. if .4 < y < .41 else 0.),
        observed_breakpoints=(.4, .41))
    # Swap the integrals. The truth endpoint corrections are below float range.
    actual = detector.expected_counts(lambda E: 1., [0, 1], 1., (0, 1),
        truth_breakpoints=[.2, .3, .32, .33, .34, .36, .4, .41, .45, .5, .6])
    assert actual[0] == pytest.approx(.01, rel=1e-7, abs=0)


@pytest.mark.parametrize('mean', [.3, .5, .67])
def test_piecewise_efficiency_partition_additivity_and_input_ownership(mean):
    boundaries = np.array([0., .4, .6, .61, 1.])
    values = np.array([.2, .4, 1., .1])

    def efficiency(y):
        assert np.isfinite(y) and 0 <= y <= 1
        return float(values[min(np.searchsorted(boundaries, y, side='right')-1, 3)])

    model = GaussianResolution(.01)
    reference = np.array([value*normal_interval((a-mean)/.01, (b-mean)/.01)
        for a, b, value in zip(boundaries[:-1], boundaries[1:], values)])
    # Duplicate and out-of-window hints are harmless and remain caller-owned.
    points = np.array([.61, .4, .6, .4, -2., 3.])
    original = points.copy()
    coarse = model.bin_probabilities(mean, [0, 1], observed_efficiency=efficiency,
                                     observed_breakpoints=points)
    fine = model.bin_probabilities(mean, boundaries, observed_efficiency=efficiency)
    np.testing.assert_allclose(fine, reference, rtol=3e-12, atol=0)
    assert coarse[0] == pytest.approx(sum(reference), rel=3e-12, abs=0)
    assert coarse[0] == pytest.approx(sum(fine), rel=3e-12, abs=0)
    np.testing.assert_array_equal(points, original)


@pytest.mark.parametrize('quenching', [ConstantQuenching(.2), LindhardQuenching(17, .11)])
@pytest.mark.parametrize('sigma', [None, 0.])
def test_declared_efficiency_boundaries_map_through_builtin_quenching(quenching, sigma):
    # The independently specified acceptance is E_nr in (.4,.41).
    lower, upper = quenching.mean_energy_keV([.4, .41])
    resolution = NoSmearing() if sigma is None else GaussianResolution(sigma)
    detector = DetectorResponse(quenching, resolution,
        observed_efficiency=CallableEfficiency(lambda y: 1. if lower < y < upper else 0.),
        observed_breakpoints=(lower, upper))
    assert detector.expected_counts(lambda E: 1., [0, 1], 1., (0, 1))[0] == pytest.approx(
        .01, rel=2e-12, abs=0)


def test_narrow_gaussian_efficiency_boundaries_map_through_nonunit_quenching():
    detector = DetectorResponse(ConstantQuenching(.2), GaussianResolution(1e-5),
        observed_efficiency=CallableEfficiency(lambda y: 1. if .08 < y < .082 else 0.),
        observed_breakpoints=(.08, .082))
    # Flat truth maps to density 1/q; Gaussian endpoint losses are negligible.
    assert detector.expected_counts(lambda E: 1., [0, 1], 1., (0, 1))[0] == pytest.approx(
        .01, rel=1e-9, abs=0)


def pipeline_limit_references(signal, background, observed):
    """Finite Poisson sum and explicit likelihood ratios, independent of backends."""
    total_signal, total_background, n = sum(signal), sum(background), int(sum(observed))

    def cdf(mean):
        term = exp(-mean)
        result = term
        for k in range(1, n+1):
            term *= mean/k
            result += term
        return result

    counting = brentq(lambda mu: cdf(total_background+mu*total_signal)/cdf(total_background)-.1,
                      0, 20, xtol=1e-13, rtol=1e-13)
    threshold = NormalDist().inv_cdf(.9)**2

    def shape(counts):
        score = lambda mu: sum(s*(ni/(b+mu*s)-1) for s, b, ni in
                               zip(signal, background, counts))
        fit = 0. if score(0) <= 0 else brentq(score, 0, 20, xtol=1e-13, rtol=1e-13)
        best = background+fit*signal

        def statistic(mu):
            means = background+mu*signal
            return 2*sum(m-b-ni*log(m/b) for m, b, ni in zip(means, best, counts))

        return brentq(lambda mu: statistic(mu)-threshold, fit, 20,
                      xtol=1e-13, rtol=1e-13)

    # For total B=2, the exact discrete median is n=2:
    assert exp(-total_background)*(1+total_background) < .5 <= cdf(total_background)
    return counting, shape(observed), shape(background)


@pytest.mark.parametrize('backend_type', [PoissonCLsCounting, BinnedPoissonProfile])
def test_discontinuous_efficiency_pipeline_limits_and_tolerance_convergence(backend_type):
    bins, signal, background, observed = [0., .405, 1.], np.array([.5, .5]), np.array([.81, 1.19]), np.array([1, 1])
    counting, shape_observed, shape_expected = pipeline_limit_references(signal, background, observed)
    detector = DetectorResponse(ConstantQuenching(1), GaussianResolution(.01),
        observed_efficiency=CallableEfficiency(lambda y: 1. if .4 < y < .41 else 0.),
        observed_breakpoints=(.4, .41))
    results = []
    for policy, roots in [(IntegrationPolicy(), RootPolicy()),
                          (IntegrationPolicy(epsrel=1e-10, limit=500), RootPolicy(xtol=1e-13, rtol=1e-13))]:
        analysis = SensitivityAnalysis(replace(detector, policy=policy), FlatBackground(.2), backend_type(roots))
        actual = [analysis.observed_limit(lambda E: 10., bins, 10., (0, 1), observed),
                  analysis.expected_limit(lambda E: 10., bins, 10., (0, 1))]
        references = [counting, counting] if backend_type is PoissonCLsCounting else [shape_observed, shape_expected]
        for result, reference in zip(actual, references):
            np.testing.assert_allclose(result.signal_counts, signal, rtol=2e-11, atol=0)
            np.testing.assert_allclose(result.background_counts, background, rtol=2e-15)
            assert result.limit.mu_upper == pytest.approx(reference, rel=2e-9, abs=0)
            if backend_type is BinnedPoissonProfile:
                assert result.limit.status == 'asymptotic_unvalidated'
                assert 'coverage unsupported/not validated' in result.limit.notes[0]
            else:
                assert result.limit.status == 'success'
                assert 'exact discrete Poisson' in result.limit.approximation
        results.append(actual)
    for default, tighter in zip(*results):
        np.testing.assert_allclose(default.signal_counts, tighter.signal_counts, rtol=2e-11, atol=0)
        assert default.limit.mu_upper == pytest.approx(tighter.limit.mu_upper, rel=2e-9, abs=0)
