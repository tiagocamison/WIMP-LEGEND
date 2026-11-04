"""Optional pre-smearing ionization acceptance for benchmark Ge detectors.

This is a phenomenological *selection* on mean visible energy, not a
microscopic electron-hole production or trigger-noise simulation. In particular,
any nonzero smooth tail below its midpoint is an explicit modeling choice.
"""
from dataclasses import dataclass
import numpy as np
from scipy.optimize import brentq
from scipy.special import ndtr
from wimp_legend._analysis import array, scalar


@dataclass(frozen=True)
class IonizationAcceptance:
    """Pre-smearing truth efficiency, evaluated in electron-equivalent keV.

    mode: 'hard' rejects energies below e50, 'erf' turns on at e50 with
      Gaussian-CDF width. Both also honor an optional strict floor.
    quenching: None for ER truth energy, or a Ge quenching model for NR truth.

    Integration points are returned in *truth* energy units for exact treatment
    of discontinuities and reliable integration of narrow turn-ons.
    """
    mode: str
    e50_keVee: float
    width_keVee: float
    floor_keVee: float = 0.0
    quenching: object = None

    def __post_init__(self):
        if type(self.mode) is not str or self.mode not in ('hard', 'erf'):
            raise ValueError('ionization acceptance must be hard or erf (none is disabled)')
        for name in ('e50_keVee', 'width_keVee', 'floor_keVee'):
            object.__setattr__(self, name, scalar(getattr(self, name), name,
                positive=name != 'floor_keVee'))
        if self.floor_keVee > self.e50_keVee:
            raise ValueError('ionization floor must not exceed e50')
        if self.quenching is not None and not callable(getattr(self.quenching, 'mean_energy_keV', None)):
            raise TypeError('quenching must implement mean_energy_keV')

    def _visible(self, truth_keV):
        truth = array(truth_keV, 'truth_keV')
        if np.any(~np.isfinite(truth)) or np.any(truth < 0):
            raise ValueError('ionization acceptance requires finite nonnegative truth energy')
        if self.quenching is None:
            return truth
        visible = array(self.quenching.mean_energy_keV(truth), 'visible_keVee')
        if visible.shape != truth.shape or np.any(~np.isfinite(visible)):
            raise ValueError('quenching must return finite visible energies with unchanged shape')
        return visible

    def __call__(self, truth_keV):
        visible = self._visible(truth_keV)
        if self.mode == 'hard':
            efficiency = (visible >= self.e50_keVee).astype(float)
        elif self.mode == 'erf':
            efficiency = ndtr((visible - self.e50_keVee) / self.width_keVee)
        else:
            raise ValueError('ionization acceptance must be hard or erf (none is disabled)')
        if self.floor_keVee > 0:
            efficiency = np.where(visible >= self.floor_keVee, efficiency, 0.)
        return float(efficiency) if efficiency.ndim == 0 else efficiency

    def integration_breakpoints_keV(self, lo, hi):
        """Truth-domain quadrature boundaries around the visible-energy turn-on."""
        lo, hi = float(lo), float(hi)
        if not (0 <= lo < hi):
            raise ValueError('invalid truth integration interval')
        visible_lo = float(self._visible(lo))
        visible_hi = float(self._visible(hi))
        if visible_hi <= visible_lo:
            raise ValueError('expected monotonic visible-energy conversion')
        candidates = [self.e50_keVee]
        if self.mode == 'erf':
            candidates += [self.e50_keVee + j*self.width_keVee
                           for j in (-8, -6, -4, -2, 2, 4, 6, 8)]
        if self.floor_keVee > 0:
            candidates.append(self.floor_keVee)
        energies = []
        for visible in sorted(set(candidates)):
            if visible_lo < visible < visible_hi:
                if self.quenching is None:
                    energies.append(float(visible))
                else:
                    energies.append(brentq(
                        lambda e: float(self._visible(e)) - visible,
                        lo, hi, xtol=1e-13, rtol=1e-12))
        return tuple(energies)
