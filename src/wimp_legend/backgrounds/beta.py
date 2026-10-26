"""Normalized truth-electron beta shapes; source BG_DATA_001 starter, hardened here.

Only *normalized, truth-electron kinetic-energy* shapes are returned, in keV^-1.
No radionuclide activities, detector effects, or physics reference spectra are supplied.

Important limitations:
- AllowedBetaShape is a reference *approximation* for allowed transitions only.
  It uses a point-Coulomb relativistic Fermi function, not a precision beta
  spectrum. Finite charge distribution, radiative, screening, exchange,
  recoil, molecular/chemical final-state and nuclear shape corrections are
  NOT implemented.
- TabulatedBetaShape is the production-relevant port for evaluated / forbidden
  spectra (e.g. BetaShape CSV). It linearly interpolates a pinned spectrum.
- For 39Ar/42Ar/42K/210Pb, do not substitute AllowedBetaShape automatically.

Point-Coulomb approximation: not a BetaShape replacement or physical certification.
"""
from __future__ import annotations

from dataclasses import dataclass
from dataclasses import field
from numbers import Real

import numpy as np
from scipy.integrate import quad
from scipy.special import loggamma, gammaln

M_E_KEV = 510.99895000
FINE_STRUCTURE_ALPHA = 7.2973525643e-3
ELECTRON_COMPTON_WAVELENGTH_FM = 386.15926764


def _energies(values):
    if any(isinstance(v, (bool, np.bool_)) for v in np.asarray(values, dtype=object).flat):
        raise ValueError("energies must exclude boolean values")
    original = np.asarray(values)
    if original.dtype.kind not in "iuf":
        raise ValueError("energies must be real numbers excluding bool")
    energies = np.asarray(original, dtype=float)
    if np.any(~np.isfinite(energies)):
        raise ValueError('Electron energies must be finite')
    return energies


def _return_shape(result, original):
    return float(result) if np.ndim(original) == 0 else result


@dataclass(frozen=True)
class AllowedBetaShape:
    """Approximate allowed beta continuum; output is a normalized dP/dT in keV^-1.

    q_keV: electron kinetic-energy endpoint (NOT atomic mass difference for beta+)
    daughter_z: nuclear charge of the daughter
    daughter_a: daughter mass number
    beta_minus: True for beta-, False for beta+
    The Fermi function uses 1.2 A^(1/3) fm as a *declared approximation*.
    """
    q_keV: float
    daughter_z: int
    daughter_a: int
    beta_minus: bool = True
    _normalization: float = field(init=False, repr=False)

    def __post_init__(self):
        if isinstance(self.q_keV, (bool, np.bool_)) or not isinstance(self.q_keV, Real) or (not np.isfinite(self.q_keV)) or self.q_keV <= 0:
            raise ValueError('Positive finite beta electron endpoint required')
        if isinstance(self.daughter_z, bool) or not isinstance(self.daughter_z, int) or not (0 <= self.daughter_z <= 118):
            raise ValueError('Invalid daughter_z')
        if isinstance(self.daughter_a, bool) or not isinstance(self.daughter_a, int) or self.daughter_a <= 0:
            raise ValueError('Invalid daughter_a')
        if not isinstance(self.beta_minus, bool):
            raise ValueError('beta_minus must be a boolean')
        if self.daughter_z * FINE_STRUCTURE_ALPHA >= 1:
            raise ValueError('Point-Coulomb gamma is undefined')
        object.__setattr__(self, '_normalization', self._compute_normalization())

    def _unnormalized(self, values):
        t = _energies(values)
        out = np.zeros_like(t, dtype=float)
        mask = (t > 0) & (t < self.q_keV)
        tt = t[mask]
        w = (tt + M_E_KEV) / M_E_KEV
        p = np.sqrt(tt * (tt + 2 * M_E_KEV)) / M_E_KEV
        z = self.daughter_z
        alpha_z = z * FINE_STRUCTURE_ALPHA
        gamma = np.sqrt(1 - alpha_z**2)
        eta = (1 if self.beta_minus else -1) * alpha_z * w / p
        r_natural = (1.2 * self.daughter_a**(1 / 3)) / ELECTRON_COMPTON_WAVELENGTH_FM
        log_f = (np.log(2 * (1 + gamma)) + 2 * (gamma - 1) * np.log(2 * p * r_natural)
                 + np.pi * eta + 2 * np.real(loggamma(gamma + 1j * eta))
                 - 2 * gammaln(2 * gamma + 1))
        out[mask] = (p * w * (self.q_keV - tt)**2 * np.exp(log_f))
        # At T=0 a beta-minus Coulomb spectrum has a finite limit: p*F
        # tends to a constant, even though p itself vanishes. A beta-plus
        # continuum and the Z=0 limit both vanish at the threshold.
        if self.beta_minus and self.daughter_z > 0:
            prefactor = (2 * (1 + gamma) * (2 * r_natural)**(2 * (gamma - 1))
                         * (2 * np.pi) * alpha_z**(2 * gamma - 1)
                         / np.exp(2 * gammaln(2 * gamma + 1)))
            out[t == 0] = prefactor * self.q_keV**2
        return _return_shape(out, values)

    def _compute_normalization(self):
        # Break at threshold, rather than relying on a midpoint quadrature there.
        result, _ = quad(lambda x: self._unnormalized(x), 0, self.q_keV,
                         epsabs=0, epsrel=2e-10, limit=250)
        if not (np.isfinite(result) and result > 0):
            raise ValueError('Beta shape is not normalizable')
        return float(result)

    def __call__(self, t_keV):
        return self._unnormalized(t_keV) / self._normalization


@dataclass(frozen=True)
class TabulatedBetaShape:
    """Frozen evaluated spectrum sampled in keV, shape density normalized per keV.

    The caller must provide the complete physical [0,Q] table. This container
    checks the zero lower edge, not the physical endpoint independently; the
    pinned Ar39 importer checks Q=565 keV. No interpolation outside the supplied
    table occurs. Raw source values remain owned tuples; source_integral records
    the full-table normalization. Keep the
    original CSV, generation tool/version, source ENSDF, branch and checksum.
    """
    energies_keV: tuple[float, ...]
    unnormalized_dndE: tuple[float, ...]
    source_id: str
    _normalization: float = field(init=False, repr=False)
    _energy_grid: np.ndarray = field(init=False, repr=False, compare=False)
    _rate_grid: np.ndarray = field(init=False, repr=False, compare=False)

    def __post_init__(self):
        energy = _energies(self.energies_keV)
        rates = _energies(self.unnormalized_dndE)
        if energy.ndim != 1 or len(energy) < 3 or rates.shape != energy.shape:
            raise ValueError('Spectrum requires at least three matched samples')
        if energy[0] != 0 or not (energy[-1] > 0) or np.any(np.diff(energy) <= 0):
            raise ValueError('Full [0,Q] domain with strictly increasing energy is required')
        if np.any(rates < 0):
            raise ValueError('Spectrum cannot be negative')
        if not isinstance(self.source_id, str) or not self.source_id.strip():
            raise ValueError('An identifying provenance record is required')
        norm = np.sum((rates[:-1]+rates[1:])*np.diff(energy)/2)
        if not np.isfinite(norm) or norm <= 0:
            raise ValueError('Spectrum is not normalizable')
        object.__setattr__(self, 'energies_keV', tuple(map(float, energy)))
        object.__setattr__(self, 'unnormalized_dndE', tuple(map(float, rates)))
        object.__setattr__(self, '_normalization', float(norm))
        # Dense-knot integration calls the scalar port repeatedly. Avoid copying
        # both source tuples per call; immutable bytes prevent write re-enabling
        # through these private NumPy views, as well as caller-owned mutation.
        object.__setattr__(self, '_energy_grid', np.frombuffer(energy.tobytes(), dtype=float))
        object.__setattr__(self, '_rate_grid', np.frombuffer(rates.tobytes(), dtype=float))

    @property
    def source_integral(self):
        return self._normalization

    @property
    def q_keV(self):
        return self.energies_keV[-1]

    def __call__(self, t_keV):
        t = _energies(t_keV)
        out = np.interp(t, self._energy_grid, self._rate_grid, left=0.0, right=0.0)
        return _return_shape(out / self._normalization, t_keV)
