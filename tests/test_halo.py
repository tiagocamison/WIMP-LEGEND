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
