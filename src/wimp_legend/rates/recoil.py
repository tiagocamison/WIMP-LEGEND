"""Generic reference elastic NREFT assembly; no nuclear data or halo assumptions.

Anand et al., PRC 89, 065501 (2014), Eqs. (40), (50), (52–53), under
NREFT_CONVENTION_ID. See docs/NREFT_CONVENTIONS.md and docs/RATE_LAYER.md.
All 32 canonical nuclear entries are required, including explicit known zeros.
"""

from dataclasses import dataclass
from math import isfinite, pi
from numbers import Real
import warnings

import numpy as np
from scipy.integrate import IntegrationWarning, quad

from ..constants import C_CM_S, GEV_MINUS2_TO_CM2, KEV_TO_GEV, SECONDS_PER_DAY
from ..interactions import (
    NREFT_CONVENTION_ID, RESPONSE_CHANNELS, WilsonCoefficients, particle_response,
)
from ..kinematics import minimum_speed_c, momentum_transfer_GeV
from ..targets import NuclearResponseDataset, ResponseKey


_EXTERNAL_Q_CHANNELS = frozenset((
    'Phi_double_prime', 'Phi_double_prime_M', 'Phi_tilde_prime',
    'Delta', 'Delta_Sigma_prime',
))
_REQUIRED_KEYS = tuple(ResponseKey(ch, t, u) for ch in RESPONSE_CHANNELS
                       for t in (0, 1) for u in (0, 1))


def _scalar(value, name, *, positive=False):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError(f'{name} must be a real scalar')
    value = float(value)
    if not isfinite(value) or (value <= 0 if positive else value < 0):
        raise ValueError(f'{name} must be finite and {"positive" if positive else "nonnegative"}')
    return value


def _array(value, name):
    result = np.asarray(value)
    if result.dtype.kind not in 'iuf':
        raise ValueError(f'{name} must contain real numbers, not bool or complex')
    result = np.asarray(result, dtype=float)
    if np.any(~np.isfinite(result)) or np.any(result < 0):
        raise ValueError(f'{name} must contain finite nonnegative values')
    return result


def _output(value):
    if np.any(~np.isfinite(value)):
        raise ValueError('rate assembly produced nonfinite values')
    return float(value) if np.ndim(value) == 0 else np.array(value, copy=True)


@dataclass(frozen=True)
class NREFTRateConfig:
    """Physics configuration: required nucleon reference mass in GeV.

    No default is defined. Changing this number changes Q and numerical
    predictions, but does not change the analytic NREFT_CONVENTION_ID.
    Store this object/value with each analysis; quadrature settings are separate.
    """

    m_N_reference_GeV: float

    def __post_init__(self):
        object.__setattr__(self, 'm_N_reference_GeV', _scalar(
            self.m_N_reference_GeV, 'm_N_reference_GeV', positive=True))


def _assembly(dataset, coefficients, config, j_chi):
    if not isinstance(dataset, NuclearResponseDataset):
        raise TypeError('dataset must be a NuclearResponseDataset')
    if dataset.metadata.convention_id != NREFT_CONVENTION_ID:
        raise ValueError('nuclear convention_id must exactly match NREFT_CONVENTION_ID')
    # Reject contradictory machine-readable declarations even when ID matches.
    if dataset.metadata.isospin_labels != ('0', '1'):
        raise ValueError("canonical nuclear isospin_labels must be ('0', '1')")
    if dataset.metadata.response_units != 'dimensionless':
        raise ValueError('canonical nuclear response_units must be dimensionless')
    if not isinstance(coefficients, WilsonCoefficients):
        raise TypeError('coefficients must be WilsonCoefficients')
    if not isinstance(config, NREFTRateConfig):
        raise TypeError('config must be NREFTRateConfig')
    spin = _scalar(j_chi, 'j_chi')
    if not (2 * spin).is_integer() or not isfinite(spin * (spin + 1)):
        raise ValueError('j_chi must have a finite spin factor and be integer or half-integer')
    for key in _REQUIRED_KEYS:
        if key not in dataset.responses:
            raise KeyError(f'{dataset.dataset_id}: missing required response {key}')
    return spin


def transition_probability(dataset, coefficients, *, q2_GeV2, V, j_chi, config):
    """Return P_tot in GeV^-4, Eq. (40), with ordered R*W and nuclear spin J.

    Q=q²/m_N_reference² is assembled here; V=v_T_perp² is dimensionless in
    units of c². q²/V broadcast. All canonical W entries must exist explicitly,
    have the exact convention ID, labels ('0','1'), and dimensionless units.
    W is evaluated through its declared source-coordinate/domain interface.
    The five composite/interference channels receive exactly one external Q.
    Signed terms/results are retained without clipping, abs(), or transposition.
    Scalar inputs return a float; arrays preserve broadcast shape, including empty.
    """
    spin = _assembly(dataset, coefficients, config, j_chi)
    q2, transverse = np.broadcast_arrays(_array(q2_GeV2, 'q2_GeV2'), _array(V, 'V'))
    with np.errstate(over='raise', invalid='raise', divide='raise'):
        Q = (np.sqrt(q2) / config.m_N_reference_GeV)**2
        if np.any((q2 > 0) & (Q == 0)):
            raise ValueError('Q underflowed; unsupported numerical scale')
        total = np.zeros_like(Q)
        for key in _REQUIRED_KEYS:
            nuclear = dataset.evaluate(channel=key.channel, tau=key.tau,
                                       tau_prime=key.tau_prime, q2_GeV2=q2)
            particle = particle_response(key.channel, coefficients, tau=key.tau,
                                         tau_prime=key.tau_prime, Q=Q, V=transverse, j_chi=spin)
            term = particle * nuclear
            if key.channel in _EXTERNAL_Q_CHANNELS:
                term = term * Q
            total = total + term
        total = total * (4*pi / (2*dataset.isotope.spin + 1))
    return _output(total)


def differential_cross_section_GeV_minus3(
    E_nr_keV, beta, *, m_chi_GeV, dataset, coefficients, j_chi, config,
):
    """Return d sigma/d E_GeV in GeV^-3, despite recoil INPUT being in keV.

    Eq. (50): m_T*P_tot/(2*pi*beta²). The supplied isotope owns m_T and J.
    E/beta arrays broadcast. Require 0 <= beta < 1 (NR physics needs beta << 1).
    Inaccessible beta < beta_min returns zero without evaluating W out of domain.
    At beta=beta_min>0, V=0 exactly; (E,beta)=(0,0) is singular and rejected.
    No relativistic-amplitude factor, composition, density or SI rescaling.
    """
    spin = _assembly(dataset, coefficients, config, j_chi)
    mass = _scalar(m_chi_GeV, 'm_chi_GeV', positive=True)
    energy, speed = np.broadcast_arrays(_array(E_nr_keV, 'E_nr_keV'), _array(beta, 'beta'))
    if np.any(speed >= 1):
        raise ValueError('beta must be below 1, in units of c')
    minimum = minimum_speed_c(energy, mass, dataset.isotope.mass_GeV)
    accessible = speed >= minimum
    if np.any(accessible & (speed == 0)):
        raise ValueError('differential cross section is undefined at E=beta=0')
    result = np.zeros_like(energy)
    if np.any(accessible):
        b, a = speed[accessible], minimum[accessible]
        with np.errstate(over='raise', invalid='raise', divide='raise'):
            V = (b - a) * (b + a)
            q = momentum_transfer_GeV(energy[accessible], dataset.isotope.mass_GeV)
            q2 = q*q
            if np.any((q > 0) & (q2 == 0)):
                raise ValueError('q2_GeV2 underflowed; unsupported numerical scale')
            probability = transition_probability(dataset, coefficients, q2_GeV2=q2,
                                                   V=V, j_chi=spin, config=config)
            result[accessible] = probability / b / b * (dataset.isotope.mass_GeV / (2*pi))
    return _output(result)


def _quadrature(beta_min, beta_max, breakpoints, epsrel, epsabs):
    lower = _scalar(beta_min, 'beta_min')
    upper = _scalar(beta_max, 'beta_max', positive=True)
    if upper >= 1:
        raise ValueError('beta_max must be below 1, in units of c')
    relative = _scalar(epsrel, 'epsrel', positive=True)
    absolute = _scalar(epsabs, 'epsabs')
    if not 50*np.finfo(float).eps < relative < 1:
        raise ValueError('epsrel must be between 50*machine_epsilon and 1')
    points = _array(tuple(breakpoints), 'breakpoints')
    if points.ndim != 1 or np.any(points <= 0) or np.any(points >= upper):
        raise ValueError('breakpoints must lie strictly inside (0, beta_max)')
    return lower, upper, np.unique(points), relative, absolute


def integrate_speed_flux(
    speed_pdf, differential_cross_section, *, beta_min, beta_max,
    breakpoints=(), epsrel=1e-8, epsabs=0.0,
):
    """Integrate f1(beta)*beta*differential_cross_section(beta) d(beta).

    Both callables accept scalar beta; cross section is in GeV^-3 and may be
    signed. Return a float in GeV^-3. f1 must be a normalized nonnegative PDF
    in d(beta), supported in [0,beta_max]; normalization is the caller's contract.
    Generic adaptive quadrature in x=beta/beta_max uses supplied breakpoints.
    No NREFT, polynomial velocity, or velocity-moment reduction is assumed.
    epsabs has the units of this integral; IntegrationWarning becomes an error.
    An empty integration interval returns zero after validating inputs.
    """
    lower, upper, points, relative, absolute = _quadrature(
        beta_min, beta_max, breakpoints, epsrel, epsabs)
    if not callable(speed_pdf) or not callable(differential_cross_section):
        raise TypeError('speed_pdf and differential_cross_section must be callable')
    if lower >= upper:
        return 0.0

    def integrand(x):
        beta = upper*x
        sample = _array(speed_pdf(beta), 'speed_pdf density')
        if sample.ndim != 0:
            raise ValueError('speed_pdf must return a scalar density')
        density = float(sample)
        cross = differential_cross_section(beta)
        if np.ndim(cross) != 0 or np.asarray(cross).dtype.kind not in 'iuf':
            raise ValueError('cross section must return a finite real scalar')
        cross = float(cross)
        if not isfinite(cross):
            raise ValueError('cross section must return a finite real scalar')
        value = (density * upper) * (beta * cross)
        if not isfinite(value):
            raise ValueError('speed-flux integrand is not representable')
        return value

    with warnings.catch_warnings(), np.errstate(over='raise', invalid='raise', divide='raise'):
        warnings.simplefilter('error', IntegrationWarning)
        integral = quad(integrand, lower/upper, 1, points=points[points > lower]/upper,
                        epsrel=relative, epsabs=absolute, limit=300)[0]
    return _output(integral)


def differential_rate_per_kg_day_keV(
    E_nr_keV, *, m_chi_GeV, rho_chi_GeV_cm3, targets_per_kg,
    dataset, coefficients, j_chi, config, speed_pdf, beta_max,
    breakpoints=(), epsrel=1e-8, epsabs=0.0,
):
    """Single-isotope reference dR/dE_nr in events/(kg day keV).

    Recoil inputs are true nuclear keV. Each array element is integrated
    independently; scalar/empty/arbitrary shape is preserved. m_chi, rho_chi
    and targets_per_kg are explicit scalars (mass positive, density/count >=0).
    No isotope fraction or mass-to-target-count conversion is implicit.
    Multiply the generic speed-flux integral by targets*rho/m, c in cm/s,
    seconds/day and (hbar*c)^2 * dE_GeV/dE_keV. beta_min>=beta_max gives zero.
    The E=0 integral is allowed if convergent; quadrature does not sample beta=0.
    epsabs refers to the natural-unit speed-flux integral, not the final rate.
    """
    spin = _assembly(dataset, coefficients, config, j_chi)
    mass = _scalar(m_chi_GeV, 'm_chi_GeV', positive=True)
    density = _scalar(rho_chi_GeV_cm3, 'rho_chi_GeV_cm3')
    count = _scalar(targets_per_kg, 'targets_per_kg')
    energy = _array(E_nr_keV, 'E_nr_keV')
    _, upper, points, relative, absolute = _quadrature(0, beta_max, breakpoints, epsrel, epsabs)
    if not callable(speed_pdf):
        raise TypeError('speed_pdf must be callable')
    minimum = minimum_speed_c(energy, mass, dataset.isotope.mass_GeV)
    result = np.zeros_like(energy)
    for index in np.ndindex(energy.shape):
        if minimum[index] >= upper:
            continue
        recoil = float(energy[index])
        def cross_section(beta):
            return differential_cross_section_GeV_minus3(
                recoil, beta, m_chi_GeV=mass, dataset=dataset,
                coefficients=coefficients, j_chi=spin, config=config)
        integral = integrate_speed_flux(speed_pdf, cross_section, beta_min=float(minimum[index]),
                                        beta_max=upper, breakpoints=points,
                                        epsrel=relative, epsabs=absolute)
        with np.errstate(over='raise', invalid='raise'):
            # Apply the area/Jacobian conversion before the large target count.
            physical_flux = integral * GEV_MINUS2_TO_CM2 * KEV_TO_GEV * C_CM_S
            result[index] = physical_flux * (density / mass) * count * SECONDS_PER_DAY
    return _output(result)
