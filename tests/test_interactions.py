"""Physics checks against the frozen Anand PRC 89, 065501 Eq. (38).

The independent reference expands bilinears into a sparse term table, without
using production coefficient accessors, formulas, channel lists, or helpers.
All numerical coefficients below are synthetic test inputs, not nuclear data.
"""

import ast
from itertools import product
from pathlib import Path

import numpy as np
import pytest

import wimp_legend.interactions as interactions
from wimp_legend.interactions import RESPONSE_CHANNELS, WilsonCoefficients, particle_response


CHANNELS = (
    "M", "Sigma_double_prime", "Sigma_prime", "Phi_double_prime",
    "Phi_double_prime_M", "Phi_tilde_prime", "Delta", "Delta_Sigma_prime",
)
ORDERED_ISOSPIN = tuple(product((0, 1), repeat=2))


def reference_eq38(raw, tau, tau_prime, q, v, spin):
    """Independent bilinear term transcription of normative Section 11.

    Each triple is (left operator, right operator, weight). Products involving
    O12/O15 are expanded, unlike the factored production implementation.
    Rational weights are directly those of Anand Eq. (38).
    """
    s = spin * (spin + 1)
    terms = {
        "M": [(1, 1, 1), (5, 5, s*q*v/3), (8, 8, s*v/3), (11, 11, s*q/3)],
        "Phi_double_prime": [
            (3, 3, q/4), (12, 12, s/12), (12, 15, -s*q/12),
            (15, 12, -s*q/12), (15, 15, s*q*q/12)],
        "Phi_double_prime_M": [(3, 1, 1), (12, 11, s/3), (15, 11, -s*q/3)],
        "Phi_tilde_prime": [(12, 12, s/12), (13, 13, s*q/12)],
        "Sigma_double_prime": [
            (10, 10, q/4), (4, 4, s/12), (4, 6, s*q/12),
            (6, 4, s*q/12), (6, 6, s*q*q/12),
            (12, 12, s*v/12), (13, 13, s*q*v/12)],
        "Sigma_prime": [
            (3, 3, q*v/8), (7, 7, v/8), (4, 4, s/12), (9, 9, s*q/12),
            (12, 12, s*v/24), (12, 15, -s*q*v/24),
            (15, 12, -s*q*v/24), (15, 15, s*q*q*v/24), (14, 14, s*q*v/24)],
        "Delta": [(5, 5, s*q/3), (8, 8, s/3)],
        "Delta_Sigma_prime": [(5, 4, s/3), (8, 9, -s/3)],
    }
    return {name: sum(raw[i][tau] * raw[j][tau_prime] * weight
                      for i, j, weight in table)
            for name, table in terms.items()}


def evaluate(channel, raw, Q=0.3, V=0.02, j_chi=0.5, tau=0, tau_prime=1):
    return particle_response(channel, WilsonCoefficients.from_isospin(raw),
                             Q=Q, V=V, j_chi=j_chi, tau=tau, tau_prime=tau_prime)


def test_channel_vocabulary():
    assert RESPONSE_CHANNELS == CHANNELS


@pytest.mark.parametrize("spin", [0, 0.5, 1])
@pytest.mark.parametrize("tau,tau_prime", ORDERED_ISOSPIN)
def test_o1_exact_baseline(spin, tau, tau_prime):
    raw = {1: (2.5, -4.)}
    q = np.array([[0], [0.2], [1.7]])
    v = np.array([0., 1e-6, 0.03, 0.5])
    for channel in CHANNELS:
        actual = evaluate(channel, raw, q, v, spin, tau, tau_prime)
        expected = raw[1][tau] * raw[1][tau_prime] if channel == "M" else 0
        np.testing.assert_array_equal(actual, np.full((3, 4), expected))


@pytest.mark.parametrize("spin", [0, 0.5, 1])
@pytest.mark.parametrize("operator", [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15])
def test_single_operator_limits(operator, spin):
    q, v, s = 0.4, 0.06, spin*(spin+1)
    weights = {
        3: {"Phi_double_prime": q/4, "Sigma_prime": q*v/8},
        4: {"Sigma_double_prime": s/12, "Sigma_prime": s/12},
        5: {"M": s*q*v/3, "Delta": s*q/3},
        6: {"Sigma_double_prime": s*q*q/12},
        7: {"Sigma_prime": v/8},
        8: {"M": s*v/3, "Delta": s/3},
        9: {"Sigma_prime": s*q/12},
        10: {"Sigma_double_prime": q/4},
        11: {"M": s*q/3},
        12: {"Phi_double_prime": s/12, "Phi_tilde_prime": s/12,
             "Sigma_double_prime": s*v/12, "Sigma_prime": s*v/24},
        13: {"Phi_tilde_prime": s*q/12, "Sigma_double_prime": s*q*v/12},
        14: {"Sigma_prime": s*q*v/24},
        15: {"Phi_double_prime": s*q*q/12, "Sigma_prime": s*q*q*v/24},
    }
    for channel in CHANNELS:
        # Opposite signs also exercise ordered off-diagonal isospin responses.
        actual = evaluate(channel, {operator: (2, -3)}, q, v, spin)
        assert actual == pytest.approx(-6 * weights[operator].get(channel, 0), abs=1e-15)


@pytest.mark.parametrize("channel,left,right,sign", [
    ("Phi_double_prime_M", 3, 1, 1),
    ("Delta_Sigma_prime", 5, 4, 1),
    ("Delta_Sigma_prime", 8, 9, -1),
])
def test_ordered_interference_signs(channel, left, right, sign):
    raw = {left: (2, 7), right: (3, 11)}
    factor = 1 if channel == "Phi_double_prime_M" else 0.25
    forward = evaluate(channel, raw, Q=0)
    reverse = evaluate(channel, raw, Q=0, tau=1, tau_prime=0)
    assert forward == sign * factor * 2 * 11
    assert reverse == sign * factor * 7 * 3
    assert forward != reverse


def test_o12_o15_cancellation_and_mixed_order():
    raw = {12: (2, 3), 15: (4, 1), 11: (5, 7)}
    # Left combination cancels at Q=1/2; the right one does not.
    assert evaluate("Phi_double_prime", raw, Q=0.5) == 0
    assert evaluate("Sigma_prime", raw, Q=0.5) == 0
    assert evaluate("Phi_double_prime_M", raw, Q=0.5) == 0
    assert evaluate("Phi_double_prime_M", raw, Q=0.5, tau=1, tau_prime=0) == 3.125
    # On crossing the cancellation, interference can be negative.
    assert evaluate("Phi_double_prime", raw, Q=1) == -0.25


def test_o4_o6_cross_terms():
    raw = {4: (2, 3), 6: (-4, 7)}
    assert evaluate("Sigma_double_prime", raw, Q=0.5) == 0
    assert evaluate("Sigma_double_prime", raw, Q=1) == -1.25


@pytest.mark.parametrize("spin", [0, 0.5, 1])
@pytest.mark.parametrize("tau,tau_prime", ORDERED_ISOSPIN)
def test_randomized_independent_eq38_all_channels(spin, tau, tau_prime):
    rng = np.random.default_rng(380065501)
    operators = (1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)
    # Includes both zero boundaries and their intersection, plus random inputs.
    q = np.concatenate(([0., 0.4], rng.uniform(0, 3, 18)))[:, None]
    v = np.concatenate(([0., 1e-6], rng.uniform(0, 0.2, 11)))[None, :]
    for _ in range(15):
        raw = {op: rng.normal(size=2) for op in operators}
        expected = reference_eq38(raw, tau, tau_prime, q, v, spin)
        for channel in CHANNELS:
            actual = evaluate(channel, raw, q, v, spin, tau, tau_prime)
            np.testing.assert_allclose(actual, np.broadcast_to(expected[channel], (20, 13)),
                                       rtol=3e-13, atol=3e-14, err_msg=channel)


@pytest.mark.parametrize("q,v,shape", [
    (0.2, 0.03, ()), (np.array(0.2), np.array(0.03), ()),
    (np.arange(4.), 0.03, (4,)), (0.2, np.arange(3.), (3,)),
    (np.arange(3.)[:, None], np.arange(4.)[None, :], (3, 4)),
    (np.array([]), 0.1, (0,)),
])
def test_broadcast_shapes_and_scalar_types(q, v, shape):
    raw = {op: (op/10, -op/7) for op in (1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)}
    for channel in CHANNELS:
        actual = evaluate(channel, raw, q, v)
        assert np.shape(actual) == shape
        assert isinstance(actual, float if shape == () else np.ndarray)
        qb, vb = np.broadcast_arrays(q, v)
        for index in np.ndindex(shape):
            point = actual if shape == () else actual[index]
            assert point == pytest.approx(evaluate(channel, raw, float(qb[index]), float(vb[index])))


def test_array_ownership():
    q, v = np.array([0., 0.4]), np.array([0., 0.06])
    q_original, v_original = q.copy(), v.copy()
    output = evaluate("M", {8: (1, 1)}, q, v)
    output[:] = 70
    np.testing.assert_array_equal(q, q_original)
    np.testing.assert_array_equal(v, v_original)
    assert evaluate("M", {8: (1, 1)}, q, v)[1] == 0.015


@pytest.mark.parametrize("field", ["Q", "V"])
@pytest.mark.parametrize("value", [-1, np.nan, np.inf, -np.inf, [0, -1], [0, np.nan],
                                        1j, [0j], True, "1", None])
def test_invalid_coordinates(field, value):
    kwargs = dict(Q=0.1, V=0.2, j_chi=0.5, tau=0, tau_prime=1)
    kwargs[field] = value
    with pytest.raises(ValueError, match=field):
        particle_response("M", WilsonCoefficients(), **kwargs)


@pytest.mark.parametrize("spin", [-0.5, 0.1, 0.50000000001, np.inf, np.nan, 1j, True, "0.5", [0.5]])
def test_invalid_spin(spin):
    with pytest.raises(ValueError, match="j_chi"):
        evaluate("M", {}, j_chi=spin)


@pytest.mark.parametrize("field", ["tau", "tau_prime"])
@pytest.mark.parametrize("value", [-1, 2, 0., True, "0"])
def test_invalid_isospin(field, value):
    with pytest.raises(ValueError, match=field):
        evaluate("M", {}, **{field: value})


def test_invalid_channel_coefficient_type_shape_and_missing_spin():
    for channel in ("unknown", "Sigma'", None, 1):
        with pytest.raises(ValueError, match="channel"):
            evaluate(channel, {})
    with pytest.raises(TypeError, match="coefficients"):
        particle_response("M", {}, Q=0, V=0, j_chi=0, tau=0, tau_prime=0)
    with pytest.raises(ValueError, match="shape mismatch"):
        evaluate("M", {}, Q=np.zeros(2), V=np.zeros(3))
    with pytest.raises(TypeError, match="j_chi"):
        particle_response("M", WilsonCoefficients(), Q=0, V=0, tau=0, tau_prime=0)


def test_numeric_range_failure_is_explicit():
    with pytest.raises((ValueError, FloatingPointError)):
        evaluate("M", {1: (1e308, 1e308)})
    with pytest.raises(FloatingPointError):
        evaluate("Phi_double_prime", {15: (1, 1)}, Q=1e308)


def test_interaction_dependency_boundary():
    # Constrain actual imports to this small layer, stdlib, and generic numerics.
    allowed = {"collections", "dataclasses", "math", "numbers", "types", "numpy",
               "wilson", "nreft"}
    directory = Path(interactions.__file__).parent
    for path in directory.glob("*.py"):
        tree = ast.parse(path.read_text())
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                assert all(alias.name.split(".")[0] in allowed for alias in node.names)
            elif isinstance(node, ast.ImportFrom):
                assert node.module.split(".")[0] in allowed
