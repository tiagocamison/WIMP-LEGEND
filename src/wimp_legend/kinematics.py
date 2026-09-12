"""Nonrelativistic two-body elastic scattering on a stationary target.

Natural units c=1: mass arguments are rest energies in GeV (m*c²), momentum
outputs are GeV (q*c), recoil energies are true nuclear-recoil keV, and speeds
are beta=v/c. All kernels broadcast NumPy inputs; scalar input returns a NumPy
scalar. No target isotope or halo parameters are read implicitly.
"""

import numpy as np

from .constants import KEV_TO_GEV


def _finite(value, name, *, positive=False):
    result = np.asarray(value, dtype=float)
    if np.any(~np.isfinite(result)) or np.any(result <= 0 if positive else result < 0):
        qualifier = "positive" if positive else "nonnegative"
        raise ValueError(f"{name} must contain finite {qualifier} values")
    return result


def reduced_mass_GeV(m_chi_GeV, m_target_GeV):
    """Reduced rest mass in GeV for finite, strictly positive masses."""
    chi = _finite(m_chi_GeV, "m_chi_GeV", positive=True)
    target = _finite(m_target_GeV, "m_target_GeV", positive=True)
    small, large = np.minimum(chi, target), np.maximum(chi, target)
    return small / (1 + small / large)


def momentum_transfer_GeV(E_nr_keV, m_target_GeV):
    """Return q*c = sqrt(2*m_target*c²*E_nr) in GeV; E_nr is in keV."""
    recoil = _finite(E_nr_keV, "E_nr_keV")
    target = _finite(m_target_GeV, "m_target_GeV", positive=True)
    return np.sqrt(2 * target * recoil * KEV_TO_GEV)


def minimum_speed_c(E_nr_keV, m_chi_GeV, m_target_GeV):
    """Return elastic v_min/c. Values beyond halo support are not clipped."""
    mu = reduced_mass_GeV(m_chi_GeV, m_target_GeV)
    return momentum_transfer_GeV(E_nr_keV, m_target_GeV) / (2 * mu)


def maximum_recoil_energy_keV(m_chi_GeV, m_target_GeV, v_max_c):
    """Return E_nr,max in keV for a supplied maximum speed in units of c.

    The caller supplies the halo endpoint. The formula assumes v_max << 1;
    v_max >= 1 is rejected to catch invalid dimensional inputs.
    """
    speed = _finite(v_max_c, "v_max_c")
    if np.any(speed >= 1):
        raise ValueError("v_max_c must be below 1 and is expressed in units of c")
    target = _finite(m_target_GeV, "m_target_GeV", positive=True)
    mu = reduced_mass_GeV(m_chi_GeV, target)
    return 2 * (mu / target) * mu * speed**2 / KEV_TO_GEV
