"""Velocity integrals used by direct-detection rate calculations."""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import interp1d


SpeedPDF = Callable[[np.ndarray], np.ndarray]


def build_velocity_integral(
    speed_pdf: SpeedPDF,
    v_max: float,
    *,
    power: float = -1.0,
    n_points: int = 10_000,
) -> Callable[[np.ndarray | float], np.ndarray]:
    r"""Return the generalized halo integral

    .. math::

        I_n(v_{\min}) = \int_{v_{\min}}^{v_{\max}} v^n f_1(v)\,dv.

    The standard mean inverse speed is obtained with ``power=-1``.
    """

    if v_max <= 0:
        raise ValueError("v_max must be positive")
    if n_points < 2:
        raise ValueError("n_points must be at least 2")

    # Avoid evaluating negative powers exactly at v=0.  For physical speed
    # distributions f_1(v) ~ v^2 near the origin, so this does not alter the
    # integral within numerical precision for the cases used here.
    v_min_grid = max(np.finfo(float).eps * v_max, 1e-15)
    v_grid = np.linspace(v_min_grid, v_max, n_points)
    integrand = np.asarray(speed_pdf(v_grid), dtype=float) * v_grid**power

    integral_grid = -cumulative_trapezoid(
        integrand[::-1],
        v_grid[::-1],
        initial=0.0,
    )[::-1]

    interpolator = interp1d(
        v_grid,
        integral_grid,
        bounds_error=False,
        fill_value=(float(integral_grid[0]), 0.0),
        assume_sorted=True,
    )

    def velocity_integral(v_threshold: np.ndarray | float) -> np.ndarray:
        values = np.asarray(v_threshold, dtype=float)
        result = np.asarray(interpolator(values), dtype=float)
        # Negative thresholds are physically equivalent to a zero threshold.
        return np.where(values <= 0.0, integral_grid[0], result)

    return velocity_integral


def mean_inverse_speed(
    speed_pdf: SpeedPDF,
    v_max: float,
    *,
    n_points: int = 10_000,
) -> Callable[[np.ndarray | float], np.ndarray]:
    """Convenience wrapper for the standard eta(v_min) halo integral."""

    return build_velocity_integral(
        speed_pdf,
        v_max,
        power=-1.0,
        n_points=n_points,
    )
