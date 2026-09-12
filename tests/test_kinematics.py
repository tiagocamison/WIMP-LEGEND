"""Elastic kinematics invariants, units, and Ge-like endpoint landmarks."""

import numpy as np
import pytest

from wimp_legend.constants import C_KM_S, KEV_TO_GEV
from wimp_legend.kinematics import (
    reduced_mass_GeV, momentum_transfer_GeV, minimum_speed_c,
    maximum_recoil_energy_keV,
)


def test_reduced_mass_symmetry_bounds_and_limits():
    rng = np.random.default_rng(237)
    first, second = 10**rng.uniform(-6, 8, (2, 500))
    mu = reduced_mass_GeV(first, second)
    np.testing.assert_array_equal(mu, reduced_mass_GeV(second, first))
    assert np.all(mu > 0)
    assert np.all(mu <= np.minimum(first, second))
    np.testing.assert_allclose(reduced_mass_GeV(first, first), first/2, rtol=1e-15)
    assert reduced_mass_GeV(1, 1e12) == pytest.approx(1, rel=2e-12)
    assert reduced_mass_GeV(1e300, 1e300) == 5e299


def test_momentum_units_and_zero_energy():
    # 50 GeV target, 10 keV recoil -> q*c = sqrt(0.001) GeV.
    assert momentum_transfer_GeV(10, 50) == pytest.approx(np.sqrt(.001), rel=1e-14)
    assert momentum_transfer_GeV(0, 50) == 0
    assert minimum_speed_c(0, 10, 50) == 0
    assert maximum_recoil_energy_keV(10, 50, 0) == 0
    assert C_KM_S == 299792.458
    assert KEV_TO_GEV == 1e-6


def test_minimum_speed_monotonic_and_sqrt_scaling():
    energy = np.geomspace(1e-6, 1e4, 300)
    for mass in [0.1, 5, 10, 40, 100, 1e6]:
        speed = minimum_speed_c(energy, mass, 76*.9315)
        assert np.all(np.diff(speed) > 0)
        np.testing.assert_allclose(minimum_speed_c(4*energy, mass, 76*.9315), 2*speed)


def test_endpoint_inversion_broadcasting():
    rng = np.random.default_rng(891)
    chi = 10**rng.uniform(-1, 4, (31, 1))
    target = 10**rng.uniform(0, 3, (1, 17))
    speed = rng.uniform(0, .005, (31, 17))
    endpoint = maximum_recoil_energy_keV(chi, target, speed)
    assert endpoint.shape == (31, 17)
    np.testing.assert_allclose(minimum_speed_c(endpoint, chi, target), speed, rtol=5e-15)
    # Maximal transfer is reversal of relative momentum in the CM frame.
    np.testing.assert_allclose(momentum_transfer_GeV(endpoint, target),
                               2*reduced_mass_GeV(chi, target)*speed, rtol=5e-15)


@pytest.mark.parametrize("mass,landmark", [(5, 4.50), (10, 15.83), (20, 50.15),
                                           (40, 134.72), (100, 354.33)])
def test_germanium_endpoint_landmarks(mass, landmark):
    actual = maximum_recoil_energy_keV(mass, 76*.9315, 810/299792.458)
    assert actual == pytest.approx(landmark, abs=.015)


def test_equal_mass_energy_transfer_and_heavy_wimp_limit():
    mass, speed = 70, .003
    assert maximum_recoil_energy_keV(mass, mass, speed) == pytest.approx(.5*mass*speed**2*1e6)
    assert maximum_recoil_energy_keV(1e15, mass, speed) == pytest.approx(2*mass*speed**2*1e6, rel=1e-12)


@pytest.mark.parametrize("invalid", [0, -1, np.nan, np.inf])
def test_invalid_masses(invalid):
    for chi, target in [(invalid, 70), (10, invalid)]:
        for call in [lambda: reduced_mass_GeV(chi, target),
                     lambda: minimum_speed_c(1, chi, target),
                     lambda: maximum_recoil_energy_keV(chi, target, .003)]:
            with pytest.raises(ValueError):
                call()
    with pytest.raises(ValueError):
        momentum_transfer_GeV(1, invalid)


@pytest.mark.parametrize("invalid", [-1, np.nan, np.inf])
def test_invalid_recoils(invalid):
    with pytest.raises(ValueError):
        momentum_transfer_GeV(invalid, 70)
    with pytest.raises(ValueError):
        minimum_speed_c(invalid, 10, 70)


@pytest.mark.parametrize("invalid", [-1, np.nan, np.inf, 1, 810])
def test_invalid_speeds(invalid):
    with pytest.raises(ValueError):
        maximum_recoil_energy_keV(10, 70, invalid)
