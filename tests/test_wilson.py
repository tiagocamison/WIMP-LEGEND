"""Analytic isospin identities and immutable coefficient input boundaries."""

from dataclasses import FrozenInstanceError

import numpy as np
import pytest

from wimp_legend.interactions import (
    CANONICAL_OPERATORS, NREFT_CONVENTION_ID, WilsonCoefficients,
)


def test_canonical_basis_and_omitted_zeros():
    assert CANONICAL_OPERATORS == (1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)
    c = WilsonCoefficients.from_isospin({1: (2.5, -3)})
    assert c.get(1, 0) == 2.5
    assert c.get(1, 1) == -3
    assert all(c.values[op] == (0, 0) for op in CANONICAL_OPERATORS if op != 1)
    assert all(pair == (0, 0) for pair in WilsonCoefficients().values.values())


@pytest.mark.parametrize("pn, expected", [
    ((6, 6), (6, 0)), ((6, 0), (3, 3)), ((0, 6), (3, -3)),
    ((-2, -2), (-2, 0)),
])
def test_analytic_half_sum_convention(pn, expected):
    c = WilsonCoefficients.from_proton_neutron({1: pn})
    assert c.values[1] == expected
    assert c.to_proton_neutron()[1] == pn


def test_multiple_operator_roundtrip_and_physical_units():
    pn = {1: (2e-7, -5e-7), 4: (-3, 7), 15: (1.25, -2.75)}
    c = WilsonCoefficients.from_proton_neutron(pn)
    for op, pair in pn.items():
        np.testing.assert_allclose(c.to_proton_neutron()[op], pair, rtol=1e-15)
    assert c.get(1, 0) == pytest.approx(-1.5e-7)
    assert c.get(1, 1) == pytest.approx(3.5e-7)


@pytest.mark.parametrize("factory", [
    WilsonCoefficients.from_isospin, WilsonCoefficients.from_proton_neutron,
])
def test_caller_input_protection_and_convention(factory):
    pair = np.array([2., 4.])
    values = {1: pair}
    c = factory(values)
    before = dict(c.values)
    pair[:] = 50
    values[4] = (8, 9)
    values.clear()
    assert dict(c.values) == before
    assert c.convention_id == NREFT_CONVENTION_ID
    with pytest.raises(TypeError):
        c.values[1] = (5, 6)
    with pytest.raises(TypeError):
        c.values[1][0] = 5
    with pytest.raises(FrozenInstanceError):
        c.values = {}
    with pytest.raises(FrozenInstanceError):
        c.convention_id = "another_convention"
    with pytest.raises(TypeError):
        c.to_proton_neutron()[1] = (7, 8)


@pytest.mark.parametrize("operator", [2, 0, 16, -1, 1., True, "1"])
def test_invalid_operator(operator):
    for factory in (WilsonCoefficients.from_isospin, WilsonCoefficients.from_proton_neutron):
        with pytest.raises(ValueError, match="operator"):
            factory({operator: (1, 0)})
    with pytest.raises(ValueError, match="operator"):
        WilsonCoefficients().get(operator, 0)


@pytest.mark.parametrize("tau", [-1, 2, 0., True, "0"])
def test_invalid_isospin(tau):
    with pytest.raises(ValueError, match="tau"):
        WilsonCoefficients().get(1, tau)


@pytest.mark.parametrize("value", [np.nan, np.inf, -np.inf, 1j, True, "1", None, [1]])
def test_invalid_coefficient(value):
    for factory in (WilsonCoefficients.from_isospin, WilsonCoefficients.from_proton_neutron):
        with pytest.raises(ValueError, match="coefficient"):
            factory({1: (0, value)})


@pytest.mark.parametrize("pair", [1, (1,), (1, 2, 3), "12", {0: 1, 1: 2}])
def test_invalid_pair(pair):
    with pytest.raises(ValueError, match="pair"):
        WilsonCoefficients.from_isospin({1: pair})


def test_non_mapping_and_numpy_scalar_inputs():
    with pytest.raises(TypeError, match="mapping"):
        WilsonCoefficients.from_isospin([(1, (0, 1))])
    c = WilsonCoefficients.from_isospin({np.int64(1): (np.float64(3), np.int64(4))})
    assert c.get(np.int64(1), np.int64(0)) == 3


def test_conversion_range_checks():
    c = WilsonCoefficients.from_proton_neutron({1: (1e308, 1e308)})
    assert c.values[1] == (1e308, 0)
    with pytest.raises(ValueError, match="finite"):
        WilsonCoefficients.from_isospin({1: (1e308, 1e308)}).to_proton_neutron()
