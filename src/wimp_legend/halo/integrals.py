r"""Distribution-independent speed moments I_n(a) = integral_a^vmax v^n f1(v) dv.

Every speed is beta=v/c and f1 is a density in d(beta). Thus moments are in
powers of beta, including eta=I_-1. To recover moments computed with km/s as
speed variable, multiply by C_KM_S**n (eta_km_s = eta_beta/C_KM_S).

The supplied PDF must be pure, nonnegative, normalized and supported inside
[0,v_max]. No SHM or interaction implementation is imported here. Scalar or
array thresholds are supported, with shape preserved and 0-d array results for
scalars. Negative thresholds are clamped to zero; +infinity returns zero, and
NaN is rejected. Integrability at zero is the caller's responsibility: for
SHM, f1(v)=O(v²), so n>-3 suffices. No small-speed cutoff regularizes a divergent
moment. For n<=-3 this API requires strictly positive thresholds, even for a
PDF that vanishes near zero. Directional rates require more than speed moments.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable
import warnings

import numpy as np
from scipy.integrate import IntegrationWarning, cumulative_trapezoid, quad

SpeedPDF = Callable[[np.ndarray | float], np.ndarray | float]


def _parameters(v_max, power, breakpoints):
    if not np.isscalar(v_max) or not np.isfinite(v_max) or not 0 < v_max < 1:
        raise ValueError("v_max must be a finite scalar in (0,1), in units of c")
    if not np.isscalar(power) or not np.isfinite(power):
        raise ValueError("moment power must be a finite real scalar")
    points = np.asarray(tuple(breakpoints), dtype=float)
    if points.ndim != 1 or np.any(~np.isfinite(points)) or np.any((points <= 0) | (points >= v_max)):
        raise ValueError("breakpoints must be finite speeds strictly inside (0,v_max)")
    return np.unique(points)


def _thresholds(v_min_c):
    thresholds = np.asarray(v_min_c, dtype=float)
    if np.any(np.isnan(thresholds)):
        raise ValueError("v_min_c must not contain NaN")
    return np.maximum(thresholds, 0)


def velocity_moment(
    speed_pdf: SpeedPDF, v_max: float, v_min_c: np.ndarray | float, *,
    n: float, breakpoints: Iterable[float] = (), epsrel: float = 1e-9,
) -> np.ndarray:
    r"""Reference adaptive quadrature for I_n(v_min_c), with explicit power n.

    Supply all known interior kinks/support boundaries through ``breakpoints``
    (e.g. halo.integration_breakpoints). Quadrature cannot reliably discover
    arbitrarily narrow features in an otherwise unknown callable.

    Integrate in x=v/v_max and factor v_max**n outside to avoid absolute-error
    tolerances swallowing small higher moments. epsrel controls relative error;
    integration warnings are raised as errors, never hidden by clipping.
    """
    points = _parameters(v_max, n, breakpoints)
    if not np.isscalar(epsrel) or not np.isfinite(epsrel) or not 50*np.finfo(float).eps < epsrel < 1:
        raise ValueError("epsrel must be finite, between 50*machine_epsilon and 1")
    thresholds = _thresholds(v_min_c)
    if n <= -3 and np.any(thresholds == 0):
        raise ValueError("n <= -3 requires a strictly positive threshold; no infrared cutoff is imposed")
    result = np.zeros_like(thresholds)

    def integrand(x):
        value = float(speed_pdf(v_max * x))
        if not np.isfinite(value) or value < 0:
            raise ValueError("speed_pdf must return finite nonnegative scalar densities")
        return x**n * value * v_max

    for index in np.ndindex(thresholds.shape):
        lower = float(thresholds[index])
        if lower >= v_max:
            continue
        with warnings.catch_warnings():
            warnings.simplefilter("error", IntegrationWarning)
            # One global error budget: a tiny boundary segment need not achieve
            # the same relative precision as the entire moment.
            value = quad(integrand, lower / v_max, 1.0,
                         points=points[points > lower] / v_max,
                         epsabs=0, epsrel=epsrel, limit=300)[0]
        result[index] = value * v_max**n
    return result


def build_velocity_integral(
    speed_pdf: SpeedPDF, v_max: float, *, power: float = -1.0,
    n_points: int = 10_000, breakpoints: Iterable[float] = (),
) -> Callable[[np.ndarray | float], np.ndarray]:
    r"""Build an approximate, nonnegative, monotone table of I_power(v_min).

    Retains the existing fast callable API. Composite trapezoids and linear
    interpolation have grid-dependent error, particularly near the endpoint;
    use ``velocity_moment`` for accuracy-controlled results. This table requires
    power >= -1 and a regular physical PDF with f1(v)=O(v²) at zero. More
    singular moments must use the reference API. No hidden state is read after
    construction: the table is a snapshot of the supplied PDF.

    The first cell is integrated adaptively, avoiding a fictitious small-speed
    cutoff. Known breakpoints are inserted into the grid. At the upper endpoint
    the interior PDF limit is used, so a hard truncation does not drop half the
    last cell. An interior PDF jump needs reference quadrature instead.
    """
    points = _parameters(v_max, power, breakpoints)
    if power < -1:
        raise ValueError("tabulation requires power >= -1; use velocity_moment for lower powers")
    if isinstance(n_points, (bool, np.bool_)) or not isinstance(n_points, (int, np.integer)) or n_points < 3:
        raise ValueError("n_points must be an integer >= 3")
    grid = np.unique(np.concatenate((np.linspace(0, v_max, n_points), points)))
    positive_grid = grid[1:].copy()
    positive_grid[-1] = np.nextafter(v_max, 0.0)
    pdf = np.asarray(speed_pdf(positive_grid), dtype=float)
    if pdf.shape != positive_grid.shape or np.any(~np.isfinite(pdf)) or np.any(pdf < 0):
        raise ValueError("speed_pdf must return finite nonnegative densities matching the input shape")
    integrand = pdf * positive_grid**power
    tails = -cumulative_trapezoid(integrand[::-1], grid[1:][::-1], initial=0)[::-1]
    # Use the reference kernel on exactly the first interval, not total minus
    # tail: subtraction of near-equal integrals would lose low-speed accuracy.
    first = float(velocity_moment(speed_pdf, grid[1], 0.0, n=power))
    table = np.concatenate(([first + tails[0]], tails))
    grid.setflags(write=False)
    table.setflags(write=False)

    def velocity_integral(v_threshold):
        values = _thresholds(v_threshold)
        return np.asarray(np.interp(values, grid, table, left=table[0], right=0.0))

    return velocity_integral


def mean_inverse_speed(
    speed_pdf: SpeedPDF, v_max: float, *, n_points: int = 10_000,
    breakpoints: Iterable[float] = (),
) -> Callable[[np.ndarray | float], np.ndarray]:
    """Existing fast eta(v_min) factory; identical to power=-1 tabulation.

    For reference quadrature use velocity_moment(..., n=-1).
    """
    return build_velocity_integral(speed_pdf, v_max, power=-1, n_points=n_points,
                                   breakpoints=breakpoints)
