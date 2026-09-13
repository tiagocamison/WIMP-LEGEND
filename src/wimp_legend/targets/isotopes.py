"""External ground-state isotope properties and explicit nuclear fractions.

No mass, spin, natural abundance, or enrichment database is bundled. Masses
are bare-nucleus rest energies m*c**2 in GeV, not atomic masses. Importers must
document any atomic-to-nuclear mass conversion in the mass source record.
"""

from dataclasses import dataclass
from math import fsum, isclose
import re
from types import MappingProxyType
from typing import Literal, Mapping

from .provenance import SourceReference, _integer, _real


@dataclass(frozen=True)
class Isotope:
    """One ground-state nuclide; identity is (Z, A), spin is J in units of hbar.

    ``sources`` must locate identity, mass_GeV, and spin individually. Structural
    validation is not a check against a periodic table or nuclear database.
    Isomers require a future explicit state identifier; do not encode them here.
    """

    symbol: str
    A: int
    Z: int
    mass_GeV: float
    spin: float
    sources: Mapping[str, SourceReference]

    def __post_init__(self):
        if not isinstance(self.symbol, str) or re.fullmatch(r"[A-Z][a-z]?", self.symbol) is None:
            raise ValueError("symbol must have element-symbol syntax")
        A, Z = _integer(self.A, "A"), _integer(self.Z, "Z")
        mass, spin = _real(self.mass_GeV, "mass_GeV"), _real(self.spin, "spin")
        if not 1 <= Z <= A:
            raise ValueError("require 1 <= Z <= A")
        if mass <= 0 or spin < 0 or not (2 * spin).is_integer():
            raise ValueError("mass must be positive; spin must be a nonnegative half-integer")
        sources = dict(self.sources)
        if set(sources) != {"identity", "mass_GeV", "spin"}:
            raise ValueError("sources must contain exactly identity, mass_GeV, spin")
        if any(not isinstance(s, SourceReference) for s in sources.values()):
            raise TypeError("isotope sources must be SourceReference records")
        for name, value in (("A", A), ("Z", Z), ("mass_GeV", mass), ("spin", spin)):
            object.__setattr__(self, name, value)
        object.__setattr__(self, "sources", MappingProxyType(sources))

    @property
    def identity(self):
        return self.Z, self.A


@dataclass(frozen=True)
class TargetComposition:
    """Ordered (isotope, fraction) pairs; iteration preserves supplied semantics.

    Fractions refer to the nuclear ensemble: number of nuclei, or mass using
    their supplied bare-nucleus masses. They are not automatically fractions
    of the measured bulk detector mass (electrons/binding require a later
    convention). Zero fractions are allowed; duplicates are always rejected.
    Sum tolerance is absolute 1e-12, relative zero; inputs are never normalized
    silently. This tolerance is numerical policy, not a physical constant.
    """

    components: tuple[tuple[Isotope, float], ...]
    fraction_kind: Literal["number", "mass"]
    source: SourceReference

    def __post_init__(self):
        pairs = self._checked_components(self.components, self.fraction_kind, self.source)
        if not isclose(fsum(f for _, f in pairs), 1.0, rel_tol=0, abs_tol=1e-12):
            raise ValueError("fractions must sum to one (absolute tolerance 1e-12)")
        object.__setattr__(self, "components", pairs)

    @staticmethod
    def _checked_components(components, fraction_kind, source):
        if fraction_kind not in ("number", "mass"):
            raise ValueError("fraction_kind must explicitly be 'number' or 'mass'")
        if not isinstance(source, SourceReference):
            raise TypeError("source must be a SourceReference")
        pairs = []
        seen = set()
        symbols, charges = {}, {}
        for isotope, fraction in components:
            if not isinstance(isotope, Isotope):
                raise TypeError("components must contain Isotope objects")
            fraction = _real(fraction, "fraction")
            if not 0 <= fraction <= 1:
                raise ValueError("fractions must lie in [0, 1]")
            if isotope.identity in seen:
                raise ValueError(f"duplicate isotope {isotope.identity}")
            seen.add(isotope.identity)
            if (symbols.setdefault(isotope.Z, isotope.symbol) != isotope.symbol
                    or charges.setdefault(isotope.symbol, isotope.Z) != isotope.Z):
                raise ValueError("inconsistent symbol/Z assignments in composition")
            pairs.append((isotope, fraction))
        if not pairs:
            raise ValueError("composition must not be empty")
        return tuple(pairs)

    @staticmethod
    def _derived_source(source, operation, components):
        return SourceReference(
            reference="wimp_legend.targets.isotopes.TargetComposition",
            locator=operation,
            kind="implementation",
            notes="Computed fractions; input_fractions retain input isotope order. "
                  "Isotope masses and their sources are retained in the composition.",
            derived_from=source,
            transformation=operation,
            input_fractions=tuple(f for _, f in components),
        )

    @classmethod
    def from_reported(cls, components, fraction_kind, source, *, normalize=False):
        """Import reported fractions; opt in with normalize=True for rounded sums.

        Default behavior is the strict constructor. Explicit normalization accepts
        any positive total, but each fraction must still lie in [0, 1]. It returns
        a derived composition retaining the reported values and source. No
        approximately normalized TargetComposition is constructed in between.
        """
        if not isinstance(normalize, bool):
            raise ValueError("normalize must be an explicit bool")
        if not normalize:
            return cls(components, fraction_kind, source)
        pairs = cls._checked_components(components, fraction_kind, source)
        total = fsum(f for _, f in pairs)
        if total == 0:
            raise ValueError("cannot normalize zero-total fractions")
        derived = cls._derived_source(source, "normalize_" + fraction_kind, pairs)
        return cls(tuple((i, f / total) for i, f in pairs), fraction_kind, derived)

    def normalized(self):
        """Return a new, explicitly normalized composition, even for an exact sum."""
        return self.from_reported(self.components, self.fraction_kind, self.source, normalize=True)

    def __iter__(self):
        return iter(self.components)

    def to_mass_fractions(self):
        """Explicitly compute w_i = x_i*m_i / sum_j(x_j*m_j)."""
        return self._convert("mass")

    def to_number_fractions(self):
        """Explicitly compute x_i = (w_i/m_i) / sum_j(w_j/m_j)."""
        return self._convert("number")

    def _convert(self, kind):
        if kind == self.fraction_kind:
            return self
        # Scaling the masses avoids overflow without altering either formula.
        scale = max(i.mass_GeV for i, _ in self) if kind == "mass" else min(i.mass_GeV for i, _ in self)
        weights = [f * (i.mass_GeV / scale if kind == "mass" else scale / i.mass_GeV) for i, f in self]
        total = fsum(weights)
        if total == 0:
            raise ValueError("fraction conversion underflowed")
        if any(f > 0 and w == 0 for (_, f), w in zip(self, weights)):
            raise ValueError("fraction conversion lost a nonzero component to underflow")
        derived = self._derived_source(self.source, self.fraction_kind + "_to_" + kind, self.components)
        return TargetComposition(tuple((i, w / total) for (i, _), w in zip(self, weights)), kind, derived)
