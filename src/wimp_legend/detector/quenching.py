"""Explicit phenomenological quenching models, not detector calibrations."""
from dataclasses import dataclass
from numbers import Integral
from typing import Protocol

import numpy as np
from .._analysis import array, output, scalar


class QuenchingModel(Protocol):
    def mean_energy_keV(self, E_nr_keV):
        """Return E_ee_keV for finite nonnegative E_nr_keV (scalar or array)."""
        ...


@dataclass(frozen=True)
class ConstantQuenching:
    q: float

    def __post_init__(self):
        q = scalar(self.q, 'q', positive=True)
        if q > 1:
            raise ValueError('q must satisfy 0 < q <= 1; zero quenching is unsupported')
        object.__setattr__(self, 'q', q)

    def mean_energy_keV(self, E_nr_keV):
        return output(self.q * array(E_nr_keV, 'E_nr_keV'))


@dataclass(frozen=True)
class LindhardQuenching:
    """Lindhard fallback with explicit Z and k; neither is a Ge default.

    epsilon=11.5 E_nr_keV Z^(-7/3), g=3 epsilon^.15+.7 epsilon^.6+epsilon,
    Q=kg/(1+kg). These numerical coefficients define the requested model,
    not empirical calibration data. See docs/DETECTOR_RESPONSE.md.
    """
    Z: int
    k: float

    def __post_init__(self):
        if isinstance(self.Z, (bool, np.bool_)) or not isinstance(self.Z, Integral) or self.Z <= 0:
            raise ValueError('Z must be a positive integer')
        object.__setattr__(self, 'Z', int(self.Z))
        object.__setattr__(self, 'k', scalar(self.k, 'k', positive=True))

    def quenching_factor(self, E_nr_keV):
        energy = array(E_nr_keV, 'E_nr_keV')
        with np.errstate(over='raise', invalid='raise', divide='raise'):
            epsilon = energy * (11.5 * self.Z**(-7/3))
            g = 3*epsilon**.15 + .7*epsilon**.6 + epsilon
            # Algebraically kg/(1+kg), avoiding overflow in kg.
            result = np.zeros_like(g)
            positive = g > 0
            result[positive] = 1 / (1 + (1 / g[positive]) / self.k)
        return output(result)

    def mean_energy_keV(self, E_nr_keV):
        energy = array(E_nr_keV, 'E_nr_keV')
        return output(energy * self.quenching_factor(energy))
