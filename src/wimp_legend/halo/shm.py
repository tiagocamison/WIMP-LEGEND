"""Standard Halo Model speed distributions.

Velocities are expressed as fractions of the speed of light.  This keeps the
halo layer independent of detector and recoil-energy unit conventions.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.special import erf


@dataclass(frozen=True)
class StandardHaloModel:
    """Truncated Maxwell-Boltzmann halo with an optional lab-frame boost.

    Parameters
    ----------
    v0:
        Maxwellian scale speed in units of c.
    vesc:
        Galactic escape speed in units of c.
    v_lab:
        Detector/lab speed relative to the Galactic frame in units of c.
    """

    v0: float
    vesc: float
    v_lab: float

    def __post_init__(self) -> None:
        if self.v0 <= 0:
            raise ValueError("v0 must be positive")
        if self.vesc <= 0:
            raise ValueError("vesc must be positive")
        if self.v_lab < 0:
            raise ValueError("v_lab must be non-negative")

    @property
    def z(self) -> float:
        return self.vesc / self.v0

    @property
    def normalization_3d(self) -> float:
        """Normalization of exp(-v^2/v0^2) inside the escape sphere."""

        z = self.z
        return (
            np.pi**1.5
            * self.v0**3
            * (erf(z) - 2.0 * z * np.exp(-z**2) / np.sqrt(np.pi))
        )

    @property
    def max_lab_speed(self) -> float:
        """Maximum speed with non-zero support in the lab frame."""

        return self.vesc + self.v_lab

    def galactic_speed_pdf(self, v: np.ndarray | float) -> np.ndarray:
        """Normalized one-dimensional speed PDF in the Galactic frame."""

        v_arr = np.asarray(v, dtype=float)
        result = np.zeros_like(v_arr)
        mask = (v_arr >= 0.0) & (v_arr < self.vesc)
        result[mask] = (
            4.0
            * np.pi
            * v_arr[mask] ** 2
            * np.exp(-(v_arr[mask] / self.v0) ** 2)
            / self.normalization_3d
        )
        return result

    def lab_speed_pdf(self, v: np.ndarray | float) -> np.ndarray:
        """Normalized speed PDF after boosting into the detector frame.

        This is the standard piecewise form for a truncated Maxwellian halo.
        For v_lab = 0, it reduces to the Galactic-frame distribution.
        """

        if self.v_lab == 0.0:
            return self.galactic_speed_pdf(v)

        v_arr = np.asarray(v, dtype=float)
        result = np.zeros_like(v_arr)

        prefactor = np.pi * self.v0**2 / (self.v_lab * self.normalization_3d)
        lower_break = max(self.vesc - self.v_lab, 0.0)
        upper_break = self.vesc + self.v_lab

        mask1 = (v_arr >= 0.0) & (v_arr < lower_break)
        result[mask1] = v_arr[mask1] * prefactor * (
            np.exp(-((v_arr[mask1] - self.v_lab) / self.v0) ** 2)
            - np.exp(-((v_arr[mask1] + self.v_lab) / self.v0) ** 2)
        )

        mask2 = (v_arr >= lower_break) & (v_arr < upper_break)
        result[mask2] = v_arr[mask2] * prefactor * (
            np.exp(-((v_arr[mask2] - self.v_lab) / self.v0) ** 2)
            - np.exp(-(self.vesc / self.v0) ** 2)
        )

        return result
