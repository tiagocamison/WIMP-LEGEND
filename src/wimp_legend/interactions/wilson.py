"""Dimensionful Wilson coefficients: Anand et al., PRC 89, 065501, Eqs. (15–18).

The normative reference is docs/NREFT_CONVENTIONS.md. No weak-scale rescaling
or nuclear inputs enter this module.
"""

from collections.abc import Mapping
from dataclasses import dataclass, field
from math import isfinite
from numbers import Integral, Real
from types import MappingProxyType

import numpy as np


NREFT_CONVENTION_ID = "anand2014_prc89_065501_half_isospin_v1"
CANONICAL_OPERATORS = (1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)


def _index(value, allowed, name):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Integral):
        raise ValueError(f"{name} must be an integer in {allowed}")
    if value not in allowed:
        raise ValueError(f"{name} must be in {allowed}")
    return int(value)


def _finite_real(value, name):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite real scalar")
    try:
        result = float(value)
    except (OverflowError, ValueError) as exc:
        raise ValueError(f"{name} must be a finite real scalar") from exc
    if not isfinite(result):
        raise ValueError(f"{name} must be a finite real scalar")
    return result


def _copy_pairs(values):
    if not isinstance(values, Mapping):
        raise TypeError("coefficients must be a mapping from operator to a pair")
    copied = {operator: (0.0, 0.0) for operator in CANONICAL_OPERATORS}
    for operator, pair in values.items():
        operator = _index(operator, CANONICAL_OPERATORS, "operator")
        if isinstance(pair, (str, bytes, Mapping)):
            raise ValueError("each operator requires a pair of real coefficients")
        try:
            first, second = pair
        except (TypeError, ValueError) as exc:
            raise ValueError("each operator requires a pair of real coefficients") from exc
        copied[operator] = (
            _finite_real(first, "coefficient"),
            _finite_real(second, "coefficient"),
        )
    return copied


@dataclass(frozen=True)
class WilsonCoefficients:
    """Immutable real c_i^tau in GeV^-2 under NREFT_CONVENTION_ID.

    ``values`` maps an operator to ``(c0, c1)``; omitted canonical operators
    are zero. Construction copies caller-owned mappings and coefficient pairs.
    Prefer the named factories to make input isospin semantics explicit.
    """

    values: Mapping[int, tuple[float, float]] = field(default_factory=dict)

    def __post_init__(self):
        object.__setattr__(self, "values", MappingProxyType(_copy_pairs(self.values)))

    @property
    def convention_id(self) -> str:
        """The fixed analytic convention; no automatic adapters are provided."""
        return NREFT_CONVENTION_ID

    @classmethod
    def from_isospin(cls, values):
        """Construct from ``{operator: (c0, c1)}``, in GeV^-2."""
        return cls(values)

    @classmethod
    def from_proton_neutron(cls, values):
        """Construct from ``{operator: (cp, cn)}``, in GeV^-2.

        c0 = (cp + cn)/2 and c1 = (cp - cn)/2, with no other scaling.
        """
        pairs = _copy_pairs(values)
        return cls({op: (cp / 2 + cn / 2, cp / 2 - cn / 2)
                    for op, (cp, cn) in pairs.items()})

    def get(self, operator: int, tau: int) -> float:
        """Return c_operator^tau in GeV^-2; reject noncanonical indices."""
        operator = _index(operator, CANONICAL_OPERATORS, "operator")
        tau = _index(tau, (0, 1), "tau")
        return self.values[operator][tau]

    def to_proton_neutron(self) -> Mapping[int, tuple[float, float]]:
        """Return immutable ``{operator: (c0+c1, c0-c1)}`` in GeV^-2.

        Values outside finite floating-point range are rejected.
        """
        return MappingProxyType(_copy_pairs({
            op: (c0 + c1, c0 - c1) for op, (c0, c1) in self.values.items()
        }))
