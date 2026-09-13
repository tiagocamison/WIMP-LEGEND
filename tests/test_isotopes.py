"""All numbers below are synthetic fixtures, not isotope database values."""

from dataclasses import FrozenInstanceError, replace

import numpy as np
import pytest

from wimp_legend.targets import Isotope, SourceReference, TargetComposition


def source():
    return SourceReference("Synthetic target tests", "tests/test_isotopes.py", "synthetic")


def isotope(A=10, mass=8.0):
    return Isotope("X", A, 4, mass, 0, {k: source() for k in ("identity", "mass_GeV", "spin")})


@pytest.mark.parametrize("field,value", [
    ("A", 0), ("A", 3), ("A", 10.0), ("A", True), ("Z", 0), ("Z", -1),
    ("Z", 11), ("Z", 4.0), ("symbol", "ge"), ("symbol", "Ge76"),
    ("mass_GeV", 0), ("mass_GeV", -1), ("mass_GeV", np.nan),
    ("mass_GeV", np.inf), ("mass_GeV", True), ("mass_GeV", [8]),
    ("spin", -0.5), ("spin", 0.3), ("spin", np.inf), ("spin", True),
])
def test_invalid_isotope(field, value):
    with pytest.raises(ValueError):
        replace(isotope(), **{field: value})


def test_isotope_properties_and_sources_are_owned():
    sources = {k: source() for k in ("identity", "mass_GeV", "spin")}
    obj = Isotope("X", np.int64(10), np.int64(4), 8, 0.5, sources)
    sources.clear()
    assert obj.identity == (4, 10)
    assert obj.spin == 0.5
    assert obj.sources["mass_GeV"].kind == "synthetic"
    with pytest.raises(TypeError):
        obj.sources["mass_GeV"] = source()
    with pytest.raises(FrozenInstanceError):
        obj.mass_GeV = 100
    with pytest.raises(ValueError):
        replace(obj, sources={})
    with pytest.raises(TypeError):
        replace(obj, sources={k: {} for k in sources.keys() | obj.sources.keys()})


@pytest.mark.parametrize("fractions", [[], [0.4, 0.5], [1.2, -0.2], [np.nan, 1], [np.inf, 0], [True, 0]])
def test_bad_compositions(fractions):
    with pytest.raises(ValueError):
        TargetComposition(tuple((isotope(10 + j), f) for j, f in enumerate(fractions)), "number", source())


def test_no_silent_normalization_or_semantics():
    # Accepted roundoff remains untouched; materially unnormalized input fails.
    fraction = 1 - 2e-13
    c = TargetComposition([(isotope(), fraction)], "number", source())
    assert c.components[0][1] == fraction
    with pytest.raises(ValueError):
        TargetComposition([(isotope(), 0.999)], "number", source())
    for kind in (None, "abundance", "percent", "NUMBER"):
        with pytest.raises(ValueError):
            TargetComposition([(isotope(), 1)], kind, source())
    with pytest.raises(TypeError):
        TargetComposition(components=[(isotope(), 1)], source=source())


def test_duplicate_detection_uses_nuclide_identity_not_properties():
    a = isotope()
    for b in (a, replace(a, mass_GeV=9), replace(a, symbol="Y")):
        with pytest.raises(ValueError, match="duplicate"):
            TargetComposition([(a, 1), (b, 0)], "mass", source())


def test_symbol_charge_consistency_and_multiple_elements():
    with pytest.raises(ValueError, match="symbol/Z"):
        TargetComposition([(isotope(), 0.5), (replace(isotope(12), symbol="Y"), 0.5)], "number", source())
    c = TargetComposition([(isotope(), 0.5), (replace(isotope(12), Z=5, symbol="Y"), 0.5)], "number", source())
    assert len(tuple(c)) == 2


def test_fraction_conversion_analytic_and_no_input_mutation():
    a, b = isotope(), isotope(12, 16)
    pairs = [[a, np.float64(0.25)], [b, 0.75]]
    c = TargetComposition(pairs, "number", source())
    pairs[0][1] = 0
    pairs.clear()
    assert list(c) == [(a, 0.25), (b, 0.75)]
    m = c.to_mass_fractions()
    assert m.fraction_kind == "mass"
    assert [f for _, f in m] == pytest.approx([1 / 7, 6 / 7])
    assert [f for _, f in m.to_number_fractions()] == pytest.approx([0.25, 0.75])
    assert m.source.derived_from is c.source
    assert m.source.transformation == "number_to_mass"
    assert m.to_mass_fractions() is m
    assert c.to_number_fractions() is c
    # The same numeric fractions mean different ensembles under different tags.
    explicit_mass = TargetComposition(c.components, "mass", source())
    assert [f for _, f in explicit_mass.to_number_fractions()] == pytest.approx([0.4, 0.6])


def test_randomized_conversion_roundtrip():
    rng = np.random.default_rng(251)  # Arbitrary deterministic test seed.
    for _ in range(30):
        fractions = rng.dirichlet(np.ones(5))
        c = TargetComposition(tuple((isotope(10 + j, m), f) for j, (m, f) in enumerate(zip(rng.uniform(1, 200, 5), fractions))), "number", source())
        result = c.to_mass_fractions().to_number_fractions()
        np.testing.assert_allclose([f for _, f in result], fractions, rtol=2e-15)


def test_pure_and_zero_fraction_components():
    c = TargetComposition([(isotope(), 0), (isotope(12), 1)], "number", source())
    assert c.to_mass_fractions().components == c.components


def test_source_immutability_and_validation():
    upstream = ["Synthetic original calculation"]
    s = replace(source(), upstream=upstream)
    upstream.clear()
    assert s.upstream == ("Synthetic original calculation",)
    for field, value in (("reference", ""), ("locator", ""), ("kind", "trusted"), ("notes", {}), ("upstream", "a string")):
        with pytest.raises(ValueError):
            replace(s, **{field: value})


@pytest.mark.parametrize("kind", ["number", "mass"])
@pytest.mark.parametrize("fractions", [(0.25, 0.7499), (0.25, 0.7501)])
def test_explicit_normalization_of_reported_fractions(kind, fractions):
    reported_source = source()
    pairs = [[isotope(), fractions[0]], [isotope(12, 16), fractions[1]]]
    original_pairs = [p.copy() for p in pairs]
    with pytest.raises(ValueError, match="sum to one"):
        TargetComposition(pairs, kind, reported_source)
    with pytest.raises(ValueError, match="sum to one"):
        TargetComposition.from_reported(pairs, kind, reported_source)
    normalized = TargetComposition.from_reported(pairs, kind, reported_source, normalize=True)
    assert pairs == original_pairs
    assert normalized.fraction_kind == kind
    assert sum(f for _, f in normalized) == pytest.approx(1, abs=2e-16)
    np.testing.assert_allclose([f for _, f in normalized], np.array(fractions) / sum(fractions))
    assert normalized.source.input_fractions == fractions
    assert normalized.source.transformation == "normalize_" + kind
    assert normalized.source.kind == "implementation"
    assert normalized.source.derived_from is reported_source
    pairs[0][1] = 0
    pairs.clear()
    assert normalized.source.input_fractions == fractions


@pytest.mark.parametrize("fraction", [1.0, 1 - 2e-13])
def test_normalized_method_preserves_original(fraction):
    original = TargetComposition([(isotope(), fraction)], "number", source())
    normalized = original.normalized()
    assert normalized is not original
    assert original.components[0][1] == fraction
    assert original.source.transformation is None
    assert normalized.components[0][1] == 1
    assert normalized.source.derived_from is original.source
    assert normalized.source.input_fractions == (fraction,)
    with pytest.raises(FrozenInstanceError):
        normalized.components = ()


@pytest.mark.parametrize("fractions", [[], [0, 0], [-0.1, 1], [1.1, 0], [np.nan, 1], [np.inf, 1], [True, 0]])
def test_normalization_rejects_invalid_inputs(fractions):
    pairs = [(isotope(10 + j), f) for j, f in enumerate(fractions)]
    with pytest.raises(ValueError):
        TargetComposition.from_reported(pairs, "number", source(), normalize=True)


def test_normalization_keeps_duplicate_and_semantic_validation():
    with pytest.raises(ValueError, match="duplicate"):
        TargetComposition.from_reported([(isotope(), 0.4), (isotope(), 0.5)], "number", source(), normalize=True)
    with pytest.raises(ValueError, match="fraction_kind"):
        TargetComposition.from_reported([(isotope(), 0.9)], "abundance", source(), normalize=True)
    with pytest.raises(ValueError, match="bool"):
        TargetComposition.from_reported([(isotope(), 0.9)], "number", source(), normalize="yes")


@pytest.mark.parametrize("kind", ["number", "mass"])
def test_provenance_chain_survives_normalization_and_both_conversions(kind):
    reported = replace(source(), upstream=("Original synthetic measurement",), notes="Rounded input")
    normalized = TargetComposition.from_reported(
        [(isotope(), 0.25), (isotope(12, 16), 0.7499)], kind, reported, normalize=True)
    converted = normalized.to_mass_fractions() if kind == "number" else normalized.to_number_fractions()
    roundtrip = converted.to_number_fractions() if kind == "number" else converted.to_mass_fractions()
    assert roundtrip.source.derived_from is converted.source
    assert converted.source.derived_from is normalized.source
    assert normalized.source.derived_from is reported
    assert reported.upstream == ("Original synthetic measurement",)
    assert reported.transformation is None
    assert converted.source.transformation == kind + "_to_" + converted.fraction_kind
    assert roundtrip.source.transformation == converted.fraction_kind + "_to_" + kind
    assert converted.source.input_fractions == tuple(f for _, f in normalized)
    assert roundtrip.source.input_fractions == tuple(f for _, f in converted)
    np.testing.assert_allclose([f for _, f in roundtrip], [f for _, f in normalized])


def test_derived_source_owns_inputs_and_rejects_incomplete_records():
    fractions = np.array([0.25, 0.7499])
    derived = SourceReference("Internal calculation", "normalization", "implementation",
                              derived_from=source(), transformation="normalize_number", input_fractions=fractions)
    fractions[:] = 0
    assert derived.input_fractions == (0.25, 0.7499)
    with pytest.raises(FrozenInstanceError):
        derived.transformation = "number_to_mass"
    for changes in ({"derived_from": None}, {"transformation": None}, {"transformation": "guess"},
                    {"kind": "primary_literature"}, {"input_fractions": [0, 0]}):
        with pytest.raises(ValueError):
            replace(derived, **changes)
