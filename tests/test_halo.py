"""Tests for halo speed distributions and velocity integrals."""

import numpy as np
from scipy.integrate import trapezoid

from wimp_legend.halo import StandardHaloModel, mean_inverse_speed


C_KM_S = 299_792.458


def benchmark_halo() -> StandardHaloModel:
    """SHM parameters matching the current O1 notebook assumptions."""

    return StandardHaloModel(
        v0=220.0 / C_KM_S,
        vesc=544.0 / C_KM_S,
        v_lab=266.0 / C_KM_S,
    )


def test_lab_speed_pdf_is_normalized() -> None:
    halo = benchmark_halo()
    v = np.linspace(0.0, halo.max_lab_speed, 100_000)
    integral = trapezoid(halo.lab_speed_pdf(v), v)
    assert np.isclose(integral, 1.0, rtol=2e-5, atol=2e-7)


def test_speed_pdf_is_nonnegative() -> None:
    halo = benchmark_halo()
    v = np.linspace(0.0, 1.2 * halo.max_lab_speed, 10_000)
    assert np.all(halo.lab_speed_pdf(v) >= 0.0)


def test_speed_pdf_vanishes_above_kinematic_support() -> None:
    halo = benchmark_halo()
    v = np.linspace(halo.max_lab_speed, 1.5 * halo.max_lab_speed, 1000)
    assert np.all(halo.lab_speed_pdf(v) == 0.0)


def test_zero_boost_reduces_to_galactic_distribution() -> None:
    halo = StandardHaloModel(v0=220.0 / C_KM_S, vesc=544.0 / C_KM_S, v_lab=0.0)
    v = np.linspace(0.0, halo.vesc, 5000)
    assert np.allclose(halo.lab_speed_pdf(v), halo.galactic_speed_pdf(v))


def test_mean_inverse_speed_is_nonnegative_and_monotonic() -> None:
    halo = benchmark_halo()
    eta = mean_inverse_speed(halo.lab_speed_pdf, halo.max_lab_speed)
    vmin = np.linspace(0.0, halo.max_lab_speed, 2000)
    values = eta(vmin)

    assert np.all(values >= -1e-12)
    assert np.all(np.diff(values) <= 1e-9)


def test_mean_inverse_speed_vanishes_above_maximum_speed() -> None:
    halo = benchmark_halo()
    eta = mean_inverse_speed(halo.lab_speed_pdf, halo.max_lab_speed)
    assert np.all(eta(np.array([halo.max_lab_speed, 1.1 * halo.max_lab_speed])) == 0.0)


# Independent checks below go beyond the original benchmark-only tests.
from dataclasses import FrozenInstanceError
from math import erf, exp, pi, sqrt

import pytest
from scipy.integrate import quad
from scipy.special import gamma, gammainc

from wimp_legend.halo import build_velocity_integral, velocity_moment


def reference_eta(halo, threshold):
    """Independent analytic eta for 0 < v_lab < vesc, in beta units."""
    a, b, z = halo.v0, halo.v_lab, halo.vesc / halo.v0
    x, y = max(threshold, 0) / a, b / a
    norm = erf(z) - 2*z*exp(-z*z)/sqrt(pi)
    if x >= z+y:
        return 0.0
    if x < z-y:
        return (erf(x+y)-erf(x-y)-4*y*exp(-z*z)/sqrt(pi))/(2*norm*b)
    return (erf(z)-erf(x-y)-2*(z+y-x)*exp(-z*z)/sqrt(pi))/(2*norm*b)


@pytest.mark.parametrize("boost", [0, 1e-12, 1, 266, 544, 700, 1500])
def test_boost_regimes(boost):
    halo = StandardHaloModel.from_km_s(v_lab_km_s=boost)
    grid = np.linspace(-halo.max_lab_speed, 1.2*halo.max_lab_speed, 10001)
    assert np.all(halo.lab_speed_pdf(grid) >= 0)
    assert np.all(halo.lab_speed_pdf(grid[grid < 0]) == 0)
    assert np.all(halo.lab_speed_pdf(grid[grid >= halo.max_lab_speed]) == 0)
    norm = quad(lambda x: float(halo.lab_speed_pdf(x * halo.max_lab_speed)) * halo.max_lab_speed,
                0, 1, points=np.asarray(halo.integration_breakpoints)/halo.max_lab_speed,
                epsabs=1e-12, epsrel=1e-12)[0]
    assert norm == pytest.approx(1, rel=2e-11)
    if boost > 544:
        assert halo.lab_speed_pdf((boost-544)/C_KM_S/2) == 0


@pytest.mark.parametrize("boost", [0, 266, 544, 700])
def test_pdf_against_independent_angular_integral(boost):
    halo = StandardHaloModel.from_km_s(v_lab_km_s=boost)
    a, b, escape = halo.v0, halo.v_lab, halo.vesc
    norm = pi**1.5 * a**3 * (erf(escape/a)-2*(escape/a)*exp(-(escape/a)**2)/sqrt(pi))
    for speed in np.linspace(0.03, 0.99, 17)*halo.max_lab_speed:
        if b == 0:
            upper = 1 if speed < escape else -1
        else:
            upper = min(1, (escape**2-speed**2-b**2)/(2*speed*b))
        expected = 0.0
        if upper > -1:
            expected = 2*pi*speed**2/norm * quad(
                lambda cosine: exp(-(speed**2+b**2+2*speed*b*cosine)/a**2), -1, upper,
                epsabs=1e-13, epsrel=1e-12)[0]
        assert float(halo.lab_speed_pdf(speed)) == pytest.approx(expected, rel=2e-12, abs=1e-12)


def test_tiny_boost_converges_without_cancellation():
    unboosted = StandardHaloModel.from_km_s(v_lab_km_s=0)
    grid = np.linspace(0, .999*unboosted.vesc, 1000)
    expected = unboosted.galactic_speed_pdf(grid)
    for boost in [1e-4, 1e-8, 1e-12]:
        halo = StandardHaloModel.from_km_s(v_lab_km_s=boost)
        np.testing.assert_allclose(halo.lab_speed_pdf(grid), expected, rtol=1e-11, atol=1e-12)


def test_small_escape_normalization():
    halo = StandardHaloModel(v0=1e-3, vesc=1e-9, v_lab=0)
    assert halo.normalization_3d == pytest.approx(4*pi*halo.vesc**3/3, rel=1e-12)
    assert velocity_moment(halo.lab_speed_pdf, halo.max_lab_speed, 0, n=0) == pytest.approx(1, rel=1e-12)


@pytest.mark.parametrize("n", [-2.9, -2, -1, 0, 1, 2, 4])
def test_zero_boost_moments_against_incomplete_gamma(n):
    halo = StandardHaloModel.from_km_s(v_lab_km_s=0)
    order = (n+3)/2
    expected = (2*pi*halo.v0**(n+3) * gamma(order) * gammainc(order, halo.z**2)
                / halo.normalization_3d)
    actual = velocity_moment(halo.lab_speed_pdf, halo.max_lab_speed, 0, n=n)
    assert actual == pytest.approx(expected, rel=2e-8)


def test_eta_against_analytic_formula_and_wrapper():
    halo = benchmark_halo()
    vmax = halo.max_lab_speed
    threshold = np.array([-1, 0, .1*vmax, halo.vesc-halo.v_lab, .7*vmax, .95*vmax, vmax, np.inf])
    exact = np.array([reference_eta(halo, v) for v in threshold])
    numerical = velocity_moment(halo.lab_speed_pdf, vmax, threshold, n=-1,
                                breakpoints=halo.integration_breakpoints)
    np.testing.assert_allclose(numerical, exact, rtol=2e-10, atol=1e-11)
    eta = mean_inverse_speed(halo.lab_speed_pdf, vmax, breakpoints=halo.integration_breakpoints)
    generic = build_velocity_integral(halo.lab_speed_pdf, vmax, power=-1,
                                      breakpoints=halo.integration_breakpoints)
    np.testing.assert_array_equal(eta(threshold), generic(threshold))
    np.testing.assert_allclose(eta(threshold), exact, rtol=3e-6, atol=1e-5)


@pytest.mark.parametrize("n", [-2.5, -1, 0, 1, 3])
def test_generic_pdf_and_moment_derivative(n):
    # Deliberately non-SHM polynomial PDF. Checks the convention and v_max**n.
    vmax = .004
    pdf = lambda v: 3*np.asarray(v)**2/vmax**3
    thresholds = np.array([[0, .2*vmax], [.6*vmax, vmax]])
    expected = 3*vmax**n/(n+3)*(1-(thresholds/vmax)**(n+3))
    actual = velocity_moment(pdf, vmax, thresholds, n=n)
    np.testing.assert_allclose(actual, expected, rtol=2e-9, atol=0)
    point, step = .4*vmax, 1e-5*vmax
    pair = velocity_moment(pdf, vmax, [point-step, point+step], n=n)
    assert (pair[1]-pair[0])/(2*step) == pytest.approx(-point**n*pdf(point), rel=1e-8)


def test_tabulation_converges_and_handles_hard_escape_endpoint():
    halo = StandardHaloModel.from_km_s(v_lab_km_s=0)
    points = np.array([0, .15, .51, .9]) * halo.max_lab_speed
    exact = velocity_moment(halo.lab_speed_pdf, halo.max_lab_speed, points, n=-1)
    errors = []
    for size in (100, 1000, 10000):
        table = mean_inverse_speed(halo.lab_speed_pdf, halo.max_lab_speed, n_points=size)
        errors.append(np.max(np.abs(table(points)-exact)))
    assert errors[1] < errors[0]/20
    assert errors[2] < errors[1]/20


def test_singular_moments_have_no_silent_cutoff():
    halo = benchmark_halo()
    with pytest.raises(ValueError, match="positive threshold"):
        velocity_moment(halo.lab_speed_pdf, halo.max_lab_speed, 0, n=-3)
    assert velocity_moment(halo.lab_speed_pdf, halo.max_lab_speed, .5*halo.max_lab_speed, n=-4) > 0
    with pytest.raises(ValueError, match="velocity_moment"):
        build_velocity_integral(halo.lab_speed_pdf, halo.max_lab_speed, power=-2)


@pytest.mark.parametrize("field,value", [("v0", 0), ("v0", np.nan), ("vesc", np.inf),
                                        ("v_lab", -1), ("v_lab", np.nan), ("v_lab", 266),
                                        ("v0", [1e-3, 2e-3])])
def test_invalid_model(field, value):
    parameters = dict(v0=220/C_KM_S, vesc=544/C_KM_S, v_lab=266/C_KM_S)
    parameters[field] = value
    with pytest.raises(ValueError):
        StandardHaloModel(**parameters)


def test_frozen_model_and_array_shapes():
    value = np.array(220/C_KM_S)
    halo = StandardHaloModel(value, 544/C_KM_S, 266/C_KM_S)
    value[...] = 1
    assert halo.v0 == 220/C_KM_S
    with pytest.raises(FrozenInstanceError):
        halo.v_lab = 0
    assert halo.lab_speed_pdf(np.zeros((2, 3))).shape == (2, 3)
    assert halo.lab_speed_pdf(0).shape == ()
    assert halo.lab_speed_pdf(np.empty((0, 2))).shape == (0, 2)
    for function in [halo.lab_speed_pdf, halo.galactic_speed_pdf]:
        with pytest.raises(ValueError):
            function(np.nan)


@pytest.mark.parametrize("kwargs", [{"n": np.nan}, {"n": 0, "epsrel": 0},
                                     {"n": 0, "breakpoints": [np.nan]}])
def test_invalid_quadrature_parameters(kwargs):
    halo = benchmark_halo()
    with pytest.raises(ValueError):
        velocity_moment(halo.lab_speed_pdf, halo.max_lab_speed, 0, **kwargs)


@pytest.mark.parametrize("size", [0, 1, 2, 2.5, True])
def test_invalid_table_size(size):
    halo = benchmark_halo()
    with pytest.raises(ValueError):
        mean_inverse_speed(halo.lab_speed_pdf, halo.max_lab_speed, n_points=size)


def test_invalid_pdf_or_threshold():
    halo = benchmark_halo()
    with pytest.raises(ValueError):
        velocity_moment(halo.lab_speed_pdf, halo.max_lab_speed, np.nan, n=0)
    with pytest.raises(ValueError):
        mean_inverse_speed(halo.lab_speed_pdf, halo.max_lab_speed)(np.nan)
    for pdf in [lambda v: -np.ones_like(v), lambda v: np.full_like(v, np.nan)]:
        with pytest.raises(ValueError):
            build_velocity_integral(pdf, .003)
        with pytest.raises(ValueError):
            velocity_moment(pdf, .003, 0, n=0)


@pytest.mark.parametrize("boost", [1e-12, 1e-8, 1e-4])
def test_reference_quadrature_tiny_boost(boost):
    halo = StandardHaloModel.from_km_s(v_lab_km_s=boost)
    assert velocity_moment(halo.lab_speed_pdf, halo.max_lab_speed, 0, n=0,
                           breakpoints=halo.integration_breakpoints) == pytest.approx(1, rel=1e-10)


@pytest.mark.parametrize("boost", [0, 266, 544, 700])
def test_second_moment_obeys_galilean_shift(boost):
    halo = StandardHaloModel.from_km_s(v_lab_km_s=boost)
    rest = StandardHaloModel.from_km_s(v_lab_km_s=0)
    rest_second = velocity_moment(rest.lab_speed_pdf, rest.max_lab_speed, 0, n=2)
    lab_second = velocity_moment(halo.lab_speed_pdf, halo.max_lab_speed, 0, n=2,
                                 breakpoints=halo.integration_breakpoints)
    assert lab_second == pytest.approx(rest_second + halo.v_lab**2, rel=2e-10)


def test_divergent_pdf_integral_is_not_silently_returned():
    from scipy.integrate import IntegrationWarning
    with pytest.raises(IntegrationWarning):
        velocity_moment(lambda v: 1/.003, .003, 0, n=-1)
