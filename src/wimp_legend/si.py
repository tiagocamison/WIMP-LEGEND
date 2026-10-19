"""Contact O1 SI normalization and a truth-spectrum adapter, outside rate engines.

Anand, Fitzpatrick & Haxton, PRC 89, 065501 (2014), Eqs. (15–18), (50),
and docs/NREFT_CONVENTIONS.md Section 17:
sigma_p[cm²] = GEV_MINUS2_TO_CM2 * mu_chi_p²/pi * (c1_p[GeV^-2])².
No weak-scale rescaling, reference nucleon mass, or nuclear data is implicit.
"""
from dataclasses import dataclass, field
from math import isfinite, pi, sqrt
from numbers import Integral, Real
from typing import Callable

import numpy as np

from ._analysis import array, scalar
from .constants import GEV_MINUS2_TO_CM2
from .interactions import WilsonCoefficients
from .kinematics import reduced_mass_GeV
from .rates import NREFTRateConfig, differential_rate_per_kg_day_keV
from .targets import NuclearResponseDataset


def _real(value, name):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real) or not isfinite(value):
        raise ValueError(f'{name} must be a finite real scalar')
    return float(value)


def _sign(value):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Integral) or value not in (-1, 1):
        raise ValueError('proton_sign must be explicitly +1 or -1')
    return int(value)


def _mu_p(m_chi_GeV, m_proton_GeV):
    chi = scalar(m_chi_GeV, 'm_chi_GeV', positive=True)
    proton = scalar(m_proton_GeV, 'm_proton_GeV', positive=True)
    return scalar(float(reduced_mass_GeV(chi, proton)), 'mu_chi_p_GeV', positive=True)


def _representable(value, nonzero, name):
    if not isfinite(value) or (nonzero and value == 0):
        raise ValueError(f'{name} is outside floating-point range')
    return float(value)


def si_proton_coefficient_GeV_minus2(sigma_p_cm2, *, m_chi_GeV, m_proton_GeV, proton_sign):
    """Invert the SI proton cross section; its global sign is explicit, not inferred.

    Real scalar inputs only. Sigma>=0, both masses>0, and proton_sign in {-1,+1}.
    Zero sigma gives zero cp; a pure-neutron model cannot be encoded by a finite
    cn/cp ratio with sigma_p=0. Use WilsonCoefficients directly for that model.
    """
    sigma = scalar(sigma_p_cm2, 'sigma_p_cm2')
    mu = _mu_p(m_chi_GeV, m_proton_GeV)
    sign = _sign(proton_sign)
    # Taking roots before dividing avoids erasing tiny but representable sigma.
    coefficient = (sqrt(pi)*sqrt(sigma)/sqrt(GEV_MINUS2_TO_CM2))/mu
    return sign * _representable(coefficient, sigma > 0, 'c1_p_GeV_minus2')


def si_proton_cross_section_cm2(c1_p_GeV_minus2, *, m_chi_GeV, m_proton_GeV):
    """Return sigma_p in cm² for a signed dimensionful proton O1 coefficient."""
    cp = _real(c1_p_GeV_minus2, 'c1_p_GeV_minus2')
    mu = _mu_p(m_chi_GeV, m_proton_GeV)
    # Do not square cp separately: that could overflow before mass factors.
    amplitude = (mu * abs(cp)) * sqrt(GEV_MINUS2_TO_CM2/pi)
    return _representable(amplitude*amplitude, cp != 0, 'sigma_p_cm2')


def si_wilson_coefficients(sigma_p_cm2, *, m_chi_GeV, m_proton_GeV,
                           fn_over_fp, proton_sign):
    """O1 only, cn=fn_over_fp*cp; use the existing half-isospin factory.

    fn_over_fp=1 is equal p/n coupling, 0 is proton-only; negative values keep
    their relative sign. There is no hidden neutron coupling or common factor.
    """
    ratio = _real(fn_over_fp, 'fn_over_fp')
    cp = si_proton_coefficient_GeV_minus2(sigma_p_cm2, m_chi_GeV=m_chi_GeV,
                                        m_proton_GeV=m_proton_GeV, proton_sign=proton_sign)
    return WilsonCoefficients.from_proton_neutron({1: (cp, ratio*cp)})


@dataclass(frozen=True, kw_only=True)
class SIO1Signal:
    """Callable reference dR/dE_nr in events/(kg day keV_nr), via generic NREFT.

    All physics parameters are explicit. targets_per_kg belongs to this single
    supplied isotope; isotope fractions are not inferred. speed_pdf is in d beta,
    beta=v/c. j_chi and config.m_N_reference_GeV retain the generic rate contract;
    m_proton_GeV is independently used for SI normalization. Callbacks must be
    pure. Combine isotope callables explicitly for a mixed-target signal.
    """
    sigma_p_ref_cm2: float
    m_chi_GeV: float
    m_proton_GeV: float
    fn_over_fp: float
    proton_sign: int
    dataset: NuclearResponseDataset
    rho_chi_GeV_cm3: float
    targets_per_kg: float
    j_chi: float
    config: NREFTRateConfig
    speed_pdf: Callable
    beta_max: float
    breakpoints: tuple[float, ...] = ()
    epsrel: float = 1e-8
    epsabs: float = 0.
    coefficients: WilsonCoefficients = field(init=False)

    def __post_init__(self):
        for name in ('sigma_p_ref_cm2', 'm_chi_GeV', 'm_proton_GeV', 'beta_max', 'epsrel'):
            object.__setattr__(self, name, scalar(getattr(self, name), name, positive=True))
        for name in ('rho_chi_GeV_cm3', 'targets_per_kg', 'j_chi', 'epsabs'):
            object.__setattr__(self, name, scalar(getattr(self, name), name))
        object.__setattr__(self, 'fn_over_fp', _real(self.fn_over_fp, 'fn_over_fp'))
        object.__setattr__(self, 'proton_sign', _sign(self.proton_sign))
        points = array(tuple(self.breakpoints), 'breakpoints')
        if points.ndim != 1:
            raise ValueError('breakpoints must be a vector')
        object.__setattr__(self, 'breakpoints', tuple(points))
        object.__setattr__(self, 'coefficients', si_wilson_coefficients(
            self.sigma_p_ref_cm2, m_chi_GeV=self.m_chi_GeV, m_proton_GeV=self.m_proton_GeV,
            fn_over_fp=self.fn_over_fp, proton_sign=self.proton_sign))
        # The public empty-array path validates all generic configuration and
        # nuclear convention/required-entry gates without sampling W or the PDF.
        differential_rate_per_kg_day_keV(np.empty(0), **self._rate_arguments())

    def _rate_arguments(self):
        return dict(m_chi_GeV=self.m_chi_GeV, rho_chi_GeV_cm3=self.rho_chi_GeV_cm3,
            targets_per_kg=self.targets_per_kg, dataset=self.dataset,
            coefficients=self.coefficients, j_chi=self.j_chi, config=self.config,
            speed_pdf=self.speed_pdf, beta_max=self.beta_max, breakpoints=self.breakpoints,
            epsrel=self.epsrel, epsabs=self.epsabs)

    @property
    def convention_id(self):
        return self.coefficients.convention_id

    def __call__(self, E_nr_keV):
        return differential_rate_per_kg_day_keV(E_nr_keV, **self._rate_arguments())

    def cross_section_limit_cm2(self, mu_upper):
        """sigma_p,upper = mu_upper*sigma_p,ref, with the backend's own CL/method.

        At CL=.90 this is sigma_p,90. Positive infinity preserves no_sensitivity;
        retain the original LimitResult to record its construction/approximation.
        This rescales cross section (not the coefficient or its square root).
        """
        if isinstance(mu_upper, (bool, np.bool_)) or not isinstance(mu_upper, Real):
            raise ValueError('mu_upper must be a nonnegative real scalar')
        mu = float(mu_upper)
        if mu == np.inf:
            return np.inf
        mu = scalar(mu, 'mu_upper')
        return _representable(mu*self.sigma_p_ref_cm2, mu > 0, 'cross section limit')
