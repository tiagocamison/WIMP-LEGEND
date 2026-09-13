"""Particle responses from Anand et al., PRC 89, 065501 (2014), Eq. (38).

See docs/NREFT_CONVENTIONS.md. Eq. (40) external Q factors belong to the
future contraction layer, not these response functions.
"""

import numpy as np

from .wilson import WilsonCoefficients, _finite_real, _index


RESPONSE_CHANNELS = (
    "M", "Sigma_double_prime", "Sigma_prime", "Phi_double_prime",
    "Phi_double_prime_M", "Phi_tilde_prime", "Delta", "Delta_Sigma_prime",
)


def _coordinate(value, name):
    array = np.asarray(value)
    if array.dtype.kind not in "iuf":
        raise ValueError(f"{name} must contain finite nonnegative real numbers")
    array = np.asarray(array, dtype=float)
    if not np.all(np.isfinite(array)) or np.any(array < 0):
        raise ValueError(f"{name} must contain finite nonnegative real numbers")
    return array


def particle_response(channel, coefficients, *, tau, tau_prime, Q, V, j_chi):
    """Evaluate ordered R_channel^(tau,tau_prime), in GeV^-4.

    Implements Anand et al., PRC 89, 065501, Eq. (38) under the
    NREFT_CONVENTION_ID exposed by ``coefficients.convention_id``.
    ``coefficients`` must be WilsonCoefficients in GeV^-2. ``Q=q^2/m_N^2``
    and ``V=v_T_perp^2`` are finite, nonnegative dimensionless inputs, with
    velocities in units of c. No numerical nucleon mass enters this kernel.
    ``j_chi`` is an explicit nonnegative integer or half-integer scalar;
    exactly integral ``2*j_chi`` is required, with no rounding tolerance.

    Q and V broadcast with NumPy rules. A zero-dimensional result is returned
    as a Python float; otherwise a new array has the broadcast shape.
    Isospin indices are ordered, including in the two mixed responses.
    The external Q factors of Eq. (40), nuclear responses, and contraction
    prefactors are NOT included. Floating-point overflow/invalid arithmetic
    raises an error; no nonfinite response is returned.
    """
    if not isinstance(channel, str) or channel not in RESPONSE_CHANNELS:
        raise ValueError(f"channel must be one of {RESPONSE_CHANNELS}")
    if not isinstance(coefficients, WilsonCoefficients):
        raise TypeError("coefficients must be WilsonCoefficients")
    tau = _index(tau, (0, 1), "tau")
    tau_prime = _index(tau_prime, (0, 1), "tau_prime")
    spin = _finite_real(j_chi, "j_chi")
    if spin < 0 or not (2 * spin).is_integer():
        raise ValueError("j_chi must be a nonnegative integer or half-integer")
    S = _finite_real(spin * (spin + 1), "j_chi*(j_chi+1)")
    Q, V = np.broadcast_arrays(_coordinate(Q, "Q"), _coordinate(V, "V"))
    a = {op: pair[tau] for op, pair in coefficients.values.items()}
    b = {op: pair[tau_prime] for op, pair in coefficients.values.items()}

    # Literal Eq. (38). In particular, do not symmetrize the mixed channels
    # or multiply the five composite/interference channels by an extra Q.
    with np.errstate(over="raise", invalid="raise"):
        if channel == "M":
            result = a[1]*b[1] + S/3 * (
                Q*V*a[5]*b[5] + V*a[8]*b[8] + Q*a[11]*b[11])
        elif channel == "Phi_double_prime":
            result = Q/4*a[3]*b[3] + S/12 * (a[12]-Q*a[15])*(b[12]-Q*b[15])
        elif channel == "Phi_double_prime_M":
            result = a[3]*b[1] + S/3 * (a[12]-Q*a[15])*b[11]
        elif channel == "Phi_tilde_prime":
            result = S/12 * (a[12]*b[12] + Q*a[13]*b[13])
        elif channel == "Sigma_double_prime":
            result = Q/4*a[10]*b[10] + S/12 * (
                a[4]*b[4] + Q*(a[4]*b[6]+a[6]*b[4]) + Q**2*a[6]*b[6]
                + V*a[12]*b[12] + Q*V*a[13]*b[13])
        elif channel == "Sigma_prime":
            result = (Q*V*a[3]*b[3] + V*a[7]*b[7])/8 + S/12 * (
                a[4]*b[4] + Q*a[9]*b[9]
                + V/2*(a[12]-Q*a[15])*(b[12]-Q*b[15])
                + Q*V/2*a[14]*b[14])
        elif channel == "Delta":
            result = S/3 * (Q*a[5]*b[5] + a[8]*b[8])
        else:  # Delta_Sigma_prime: ordered positive c5*c4, negative c8*c9.
            result = S/3 * (a[5]*b[4] - a[8]*b[9])
    if not np.all(np.isfinite(result)):
        raise ValueError("particle response is outside finite floating-point range")
    result = np.broadcast_to(result, Q.shape)
    return float(result) if result.ndim == 0 else result.copy()
