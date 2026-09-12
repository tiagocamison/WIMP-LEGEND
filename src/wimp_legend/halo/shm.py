"""Isotropic Maxwellian truncated in the Galactic frame, then Galilean boosted.

All speeds and PDF arguments are dimensionless beta = speed/c. The PDFs are
probability densities with respect to d(beta), not d(speed in km/s).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.special import exprel, gammainc

from ..constants import C_KM_S


@dataclass(frozen=True)
class StandardHaloModel:
    r"""Physical configuration for f_G(u) proportional to exp(-u²/v0²).

    ``v0``, ``vesc``, and ``v_lab`` are scalar speeds in units of c. These
    existing constructor names are retained for compatibility; use
    ``from_km_s`` at a dimensional input boundary. v0 is the Maxwellian scale,
    not the one-component dispersion (which is v0/sqrt(2)). Truncation is a
    hard cut at Galactic speed vesc, not a cut imposed after boosting.

    No density, annual modulation, focusing, or interaction parameters are
    included. A fixed v_lab represents a snapshot, not an annual average.
    The Galilean model is physically applicable only at nonrelativistic speeds.
    """

    v0: float
    vesc: float
    v_lab: float

    def __post_init__(self) -> None:
        for name in ("v0", "vesc", "v_lab"):
            value = np.asarray(getattr(self, name), dtype=float)
            if value.ndim != 0 or not np.isfinite(value):
                raise ValueError(f"{name} must be a finite scalar speed in units of c")
            # Copy scalar inputs so even a caller-owned 0-d array cannot mutate us.
            object.__setattr__(self, name, float(value))
        if not 0 < self.v0 < 1 or self.vesc <= 0 or self.v_lab < 0:
            raise ValueError("require 0 < v0 < 1, vesc > 0, and v_lab >= 0 (units of c)")
        if self.max_lab_speed >= 1:
            raise ValueError("vesc + v_lab must be below c; use from_km_s for km/s inputs")
        if not np.isfinite(self.normalization_3d) or self.normalization_3d <= 0:
            raise ValueError("halo normalization is not representable in float64")

    @classmethod
    def from_km_s(
        cls, *, v0_km_s: float = 220.0, vesc_km_s: float = 544.0,
        v_lab_km_s: float = 266.0,
    ) -> StandardHaloModel:
        """Construct from explicitly dimensional speeds; defaults are benchmarks."""
        return cls(v0_km_s / C_KM_S, vesc_km_s / C_KM_S, v_lab_km_s / C_KM_S)

    @property
    def z(self) -> float:
        return self.vesc / self.v0

    @property
    def normalization_3d(self) -> float:
        """Integral of exp(-u²/v0²) over the escape sphere in d³(beta)."""
        # P(3/2,z²) = erf(z) - 2*z*exp(-z²)/sqrt(pi), without cancellation
        # between two nearly equal terms for a small escape-to-scale ratio.
        return float(np.pi**1.5 * self.v0**3 * gammainc(1.5, self.z**2))

    @property
    def max_lab_speed(self) -> float:
        """Upper support endpoint in units of c."""
        return self.vesc + self.v_lab

    @property
    def integration_breakpoints(self) -> tuple[float, ...]:
        """Interior support boundary/kink in units of c, for quadrature."""
        point = abs(self.vesc - self.v_lab)
        # A kink effectively coincident with an endpoint makes QUADPACK split
        # into an unresolvable interval. Leave such a kink to global adaptive
        # quadrature, instead of requesting a machine-precision-width segment.
        margin = 64 * np.finfo(float).eps * self.max_lab_speed
        return (point,) if margin < point < self.max_lab_speed - margin else ()

    def galactic_speed_pdf(self, v: np.ndarray | float) -> np.ndarray:
        """Density with respect to d(beta); scalar input gives a 0-d array.

        Negative speeds and speeds outside support give zero. NaN is rejected;
        either infinity is outside support. The value at the escape cut is zero.
        """
        values = np.asarray(v, dtype=float)
        if np.any(np.isnan(values)):
            raise ValueError("speed must not contain NaN")
        result = np.zeros_like(values)
        mask = (values >= 0) & (values < self.vesc)
        speed = values[mask]
        result[mask] = (4 * np.pi * speed**2 / self.normalization_3d
                        * np.exp(-(speed / self.v0)**2))
        return result

    def lab_speed_pdf(self, v: np.ndarray | float) -> np.ndarray:
        """Boosted speed density in d(beta), with the same input policy as above.

        Angular integration runs over |v-v_lab| < u < min(v+v_lab, vesc).
        This support condition also covers v_lab >= vesc. exprel/expm1 remove
        subtraction cancellation at tiny boosts and near the upper endpoint.
        """
        if self.v_lab == 0:
            return self.galactic_speed_pdf(v)
        values = np.asarray(v, dtype=float)
        if np.any(np.isnan(values)):
            raise ValueError("speed must not contain NaN")
        result = np.zeros_like(values)
        full = (values >= 0) & (values < self.vesc - self.v_lab)
        speed = values[full]
        delta = 4 * (speed / self.v0) * (self.v_lab / self.v0)
        result[full] = (4 * np.pi * speed**2 / self.normalization_3d
                        * np.exp(-((speed - self.v_lab) / self.v0)**2)
                        * exprel(-delta))

        partial = ((values >= abs(self.vesc - self.v_lab))
                   & (values < self.max_lab_speed) & (values > 0))
        speed = values[partial]
        lower = np.abs(speed - self.v_lab)
        delta = ((self.vesc - lower) / self.v0) * ((self.vesc + lower) / self.v0)
        result[partial] = (np.pi * self.v0**2 * speed
                           / (self.v_lab * self.normalization_3d)
                           * np.exp(-(lower / self.v0)**2) * (-np.expm1(-delta)))
        return result
