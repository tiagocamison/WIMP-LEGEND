"""Immutable source records, not certificates of physical correctness."""

from dataclasses import dataclass
from math import isfinite
from numbers import Integral, Real
from typing import Literal


def _text(value, name):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be nonempty text")
    return value


def _real(value, name):
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a real scalar")
    value = float(value)
    if not isfinite(value):
        raise ValueError(f"{name} must be finite")
    return value


def _integer(value, name):
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise ValueError(f"{name} must be an integer")
    return int(value)


@dataclass(frozen=True)
class SourceReference:
    """Identify the immediate source and its upstream calculation separately.

    ``reference`` should include publication/version or repository/commit;
    ``locator`` identifies the table, equation, file, rows, or coefficient order.
    ``upstream`` records original calculations as citations, when traced.
    Composition derivations use a single ``derived_from`` link to the intact
    previous source, a ``transformation`` tag, and owned ``input_fractions`` in
    composition isotope order. Tags also identify the input fraction semantics.
    No transformation means directly supplied values, not necessarily measured
    or validated ones. This is a linear audit trail, not a provenance graph.
    A synthetic or unresolved source must say so; no status implies validation.
    See docs/TARGETS_NUCLEAR_RESPONSES.md for the audited reference chain.
    """

    reference: str
    locator: str
    kind: Literal["primary_literature", "implementation", "synthetic", "unresolved"]
    upstream: tuple[str, ...] = ()
    notes: str = ""
    derived_from: "SourceReference | None" = None
    transformation: Literal["normalize_number", "normalize_mass", "number_to_mass", "mass_to_number"] | None = None
    input_fractions: tuple[float, ...] = ()

    def __post_init__(self):
        _text(self.reference, "reference")
        _text(self.locator, "locator")
        if self.kind not in ("primary_literature", "implementation", "synthetic", "unresolved"):
            raise ValueError("unknown source kind")
        if isinstance(self.upstream, str):
            raise ValueError("upstream must be a sequence of citations")
        object.__setattr__(self, "upstream", tuple(_text(s, "upstream citation") for s in self.upstream))
        if not isinstance(self.notes, str):
            raise ValueError("notes must be text")
        fractions = tuple(_real(v, "input fraction") for v in self.input_fractions)
        if self.transformation is None:
            if self.derived_from is not None or fractions:
                raise ValueError("derived provenance requires a transformation")
        else:
            if self.transformation not in ("normalize_number", "normalize_mass", "number_to_mass", "mass_to_number"):
                raise ValueError("unknown composition transformation")
            if self.kind != "implementation" or not isinstance(self.derived_from, SourceReference):
                raise ValueError("derived values require implementation kind and an upstream source")
            if not fractions or any(v < 0 or v > 1 for v in fractions) or not any(fractions):
                raise ValueError("derived provenance requires nonnegative, nonzero input fractions in [0, 1]")
        object.__setattr__(self, "input_fractions", fractions)
