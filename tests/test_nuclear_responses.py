"""Synthetic numerical tests only; no nuclear-structure data are validated."""

from dataclasses import FrozenInstanceError, replace

import numpy as np
import pytest

from wimp_legend.targets import (
    Isotope, TargetComposition, SourceReference, MomentumVariable, ResponseKey,
    ResponseMetadata, PolynomialResponse, TabulatedResponse, NuclearResponseDataset,
)


def source():
    return SourceReference("Synthetic response tests", "tests/test_nuclear_responses.py; all numerical fixtures", "synthetic")


def dataset(momentum=None):
    isotope = Isotope("X", 10, 4, 8, 0, {k: source() for k in ("identity", "mass_GeV", "spin")})
    metadata = ResponseMetadata(
        convention_id="synthetic-convention-v1",
        basis="synthetic-bilinear-v1", isospin_labels=("0", "1"),
        isospin_definition="Synthetic labels only; not a physical isospin convention",
        normalization="Raw synthetic polynomial; no spin average, 4*pi, q/m_N or interference factors",
        response_units="dimensionless", momentum=momentum or MomentumVariable("q2_GeV2"), source=source(),
    )
    return NuclearResponseDataset("synthetic-X10", isotope, metadata, {
        ResponseKey("M", 0, 1): PolynomialResponse([1, -2, 3], 2, [0, 4], source()),
        ResponseKey("Delta_Sigma_prime", 1, 0): TabulatedResponse([0, 1, 4], [-1, -3, -9], source()),
    })


def evaluate(d, q2, **kwargs):
    return d.evaluate(channel="M", tau=0, tau_prime=1, q2_GeV2=q2, **kwargs)


def test_polynomial_scalar_vector_and_matrix():
    d = dataset()
    assert isinstance(evaluate(d, 0), float)
    assert evaluate(d, 0) == 1
    for x in (0.25, [0, 0.25, 1, 4], np.array([[0.1, 0.3], [1, 2]]), np.empty((0, 2))):
        a = np.asarray(x)
        result = evaluate(d, x)
        np.testing.assert_allclose(result, (1 - 2 * a + 3 * a**2) * np.exp(-2 * a), rtol=1e-14)
        assert np.shape(result) == a.shape


def test_table_signed_values_and_boundaries():
    d = dataset()
    x = np.array([[0, 0.25], [1, 4]])
    result = d.evaluate(channel="Delta_Sigma_prime", tau=1, tau_prime=0, q2_GeV2=x)
    np.testing.assert_allclose(result, -1 - 2 * x)
    assert d.evaluate(channel="Delta_Sigma_prime", tau=1, tau_prime=0, q2_GeV2=0) == -1
    result[:] = 100
    assert d.responses[ResponseKey("Delta_Sigma_prime", 1, 0)].values[0] == -1


@pytest.mark.parametrize("variable,expected", [("q2_GeV2", 0.25), ("q_GeV", 0.5), ("y", 1), ("u", 2)])
def test_momentum_definitions(variable, expected):
    # b=4 GeV^-1 is deliberately synthetic. y=q²b²/4, u=q²b²/2.
    m = MomentumVariable(variable, 4, source()) if variable in ("y", "u") else MomentumVariable(variable)
    assert m.from_q2(0.25) == expected
    np.testing.assert_allclose(evaluate(dataset(m), 0.25), (1 - 2 * expected + 3 * expected**2) * np.exp(-2 * expected))


def test_equivalent_q2_and_y_parametrizations():
    # Independent analytic reparameterization with b=6, hence y=9*q².
    d_y = dataset(MomentumVariable("y", 6, source()))
    d_q = dataset()
    p = PolynomialResponse([1, -18, 243], 18, [0, 4 / 9], source())
    d_q = replace(d_q, responses={ResponseKey("M", 0, 1): p})
    q2 = np.linspace(0, 0.4, 31)
    np.testing.assert_allclose(evaluate(d_y, q2), evaluate(d_q, q2), rtol=2e-14)


@pytest.mark.parametrize("value", [-1, np.nan, np.inf, -np.inf, [0, -1], 1j, True, "0.1"])
def test_bad_q2(value):
    with pytest.raises(ValueError):
        evaluate(dataset(), value)


@pytest.mark.parametrize("indices", [(-1, 0), (0, 2), (True, 0), (0, 1.0), ("p", 0), (np.nan, 0)])
def test_invalid_isospin(indices):
    with pytest.raises(ValueError):
        dataset().evaluate(channel="M", tau=indices[0], tau_prime=indices[1], q2_GeV2=0)


def test_missing_is_not_zero_or_transpose():
    d = dataset()
    for channel, tau, tau_prime in [("Sigma_prime", 0, 0), ("M", 1, 0), ("Sigma_prime_Delta", 1, 0)]:
        with pytest.raises(KeyError, match="missing response"):
            d.evaluate(channel=channel, tau=tau, tau_prime=tau_prime, q2_GeV2=0)
    empty = replace(d, responses={})
    with pytest.raises(KeyError):
        evaluate(empty, 0)
    zero = replace(d, responses={ResponseKey("M", 0, 1): PolynomialResponse([0], 0, [0, 1], source())})
    assert evaluate(zero, 0.5) == 0


@pytest.mark.parametrize("response", [PolynomialResponse([1], 0, [0.1, 1], source()), TabulatedResponse([0.1, 1], [1, 2], source())])
def test_no_extrapolation(response):
    for value in (0, 1.01):
        with pytest.raises(ValueError, match="outside"):
            response.evaluate(value)


def test_conventions_and_provenance_preserved_without_conversion():
    d = dataset(MomentumVariable("y", 4, source()))
    meta = replace(d.metadata, basis="alternative-basis", isospin_labels=["p", "n"],
                   isospin_definition="0=p, 1=n; explicit synthetic projectors",
                   response_units="synthetic-units", normalization="Already contains a synthetic factor; return unchanged")
    other = replace(d, metadata=meta, dataset_id="alternative")
    assert other.metadata == meta
    assert other.metadata.momentum.b_source == source()
    assert other.responses[ResponseKey("M", 0, 1)].source == source()
    assert evaluate(other, 0.25) == evaluate(d, 0.25)


def test_no_caller_owned_mutable_storage():
    coeffs = np.array([1., 2.])
    domain = np.array([0., 4.])
    nodes, values = np.array([0., 4.]), np.array([1., 9.])
    p = PolynomialResponse(coeffs, 0, domain, source())
    t = TabulatedResponse(nodes, values, source())
    mapping = {ResponseKey("M", 0, 1): p}
    d = replace(dataset(), responses=mapping)
    coeffs[:] = domain[:] = nodes[:] = values[:] = 99
    mapping.clear()
    assert evaluate(d, 1) == t.evaluate(1) == 3
    with pytest.raises(TypeError):
        d.responses[ResponseKey("M", 0, 1)] = t
    with pytest.raises(FrozenInstanceError):
        p.coefficients = (7,)
    labels = ["0", "1"]
    meta = replace(d.metadata, isospin_labels=labels)
    labels[0] = "p"
    assert meta.isospin_labels == ("0", "1")


def test_composition_independence_and_multiple_datasets():
    d = dataset()
    original = evaluate(d, [0, 1])
    other_isotope = replace(d.isotope, A=12, mass_GeV=10)
    for x in (0, 0.25, 1):
        c = TargetComposition([(d.isotope, x), (other_isotope, 1 - x)], "number", source())
        c.to_mass_fractions()
        np.testing.assert_array_equal(evaluate(d, [0, 1]), original)
    alt = replace(d, dataset_id="second-calculation", responses={ResponseKey("M", 0, 1): PolynomialResponse([7], 0, [0, 4], source())})
    assert alt.isotope is d.isotope
    assert evaluate(alt, 0) == 7 and evaluate(d, 0) == 1


@pytest.mark.parametrize("coefficients,decay,domain", [([], 0, [0, 1]), ([[1]], 0, [0, 1]), ([np.nan], 0, [0, 1]), ([1j], 0, [0, 1]), ([1], -1, [0, 1]), ([1], np.inf, [0, 1]), ([1], 0, [1, 1]), ([1], 0, [-1, 1]), ([1], 0, [0, np.inf])])
def test_invalid_polynomial(coefficients, decay, domain):
    with pytest.raises(ValueError):
        PolynomialResponse(coefficients, decay, domain, source())


@pytest.mark.parametrize("nodes,values", [([0], [1]), ([0, 0], [1, 2]), ([1, 0], [1, 2]), ([-1, 0], [1, 2]), ([0, 1], [1]), ([0, 1], [1, np.inf])])
def test_invalid_tables(nodes, values):
    with pytest.raises(ValueError):
        TabulatedResponse(nodes, values, source())


def test_invalid_metadata_and_nonfinite_output():
    for args in [("other",), ("y",), ("u", -1, source()), ("y", 1, None), ("q2_GeV2", 1, source())]:
        with pytest.raises(ValueError):
            MomentumVariable(*args)
    for field, value in [("basis", ""), ("normalization", ""), ("isospin_definition", ""), ("response_units", ""), ("isospin_labels", ["0", "0"]), ("isospin_labels", "01")]:
        with pytest.raises(ValueError):
            replace(dataset().metadata, **{field: value})
    with pytest.raises(TypeError):
        replace(dataset(), responses={ResponseKey("M", 0, 0): lambda x: x})
    with pytest.raises(ValueError, match="nonfinite"):
        PolynomialResponse([1e308, 1e308], 0, [0, 4], source()).evaluate(4)


@pytest.mark.parametrize("domain", [None, (0, None)])
def test_analytic_response_without_upper_bound(domain):
    p = PolynomialResponse([1, 2], 0, domain, source())
    assert p.domain == (0, None)
    assert p.validated_range is None
    assert p.evaluate(1e10) == 1 + 2e10
    d = replace(dataset(), responses={ResponseKey("M", 0, 1): p})
    assert evaluate(d, 1e10) == 1 + 2e10
    for invalid in (-1e-20, np.inf, np.nan):
        with pytest.raises(ValueError):
            p.evaluate(invalid)


def test_unbounded_domain_can_keep_a_nonzero_lower_bound():
    p = PolynomialResponse([1], 0, (0.1, None), source())
    assert p.evaluate(1e10) == 1
    with pytest.raises(ValueError, match="domain"):
        p.evaluate(0)


@pytest.mark.parametrize("tabulated", [False, True])
def test_source_supported_range_is_separate_owned_and_enforced(tabulated):
    supported = np.array([0.25, 0.75])
    evidence = replace(source(), locator="Synthetic range-validation record")
    if tabulated:
        response = TabulatedResponse([0, 1], [1, 3], source(), supported, evidence)
    else:
        response = PolynomialResponse([1, 2], 0, None, source(), supported, evidence)
    supported[:] = 99
    assert response.validated_range == (0.25, 0.75)
    assert response.validated_range_source is evidence
    assert response.domain == ((0, 1) if tabulated else (0, None))
    np.testing.assert_allclose(response.evaluate([0.25, 0.75]), [1.5, 2.5])
    with pytest.raises(FrozenInstanceError):
        response.validated_range = None
    for invalid in (0.1, 0.9):
        with pytest.raises(ValueError, match="validated_range"):
            response.evaluate(invalid)
    d = replace(dataset(), responses={ResponseKey("M", 0, 1): response})
    with pytest.raises(ValueError, match="validated_range"):
        evaluate(d, 0.9)


@pytest.mark.parametrize("bounds", [(-1, 1), (0.5, 0.5), (0, np.inf), (0, 3), (0, np.nan)])
def test_invalid_validated_ranges(bounds):
    for response in (PolynomialResponse([1], 0, (0, 2), source()),
                     TabulatedResponse([0, 2], [1, 1], source())):
        with pytest.raises(ValueError):
            replace(response, validated_range=bounds, validated_range_source=source())


def test_validated_range_requires_source_and_source_requires_range():
    p = PolynomialResponse([1], 0, None, source())
    with pytest.raises(TypeError, match="source"):
        replace(p, validated_range=(0, 1))
    with pytest.raises(ValueError, match="requires"):
        replace(p, validated_range_source=source())


def test_table_stays_finite_and_rejects_negative_coordinates():
    table = TabulatedResponse([0, 1], [1, 2], source())
    for x in (-1e-20, 1.01, 1e10):
        with pytest.raises(ValueError, match="domain"):
            table.evaluate(x)
    with pytest.raises(ValueError):
        TabulatedResponse([0, np.inf], [1, 2], source())


def test_convention_identity_is_preserved_distinct_and_never_converted():
    d = dataset()
    alternate = replace(d, metadata=replace(d.metadata, convention_id="synthetic-convention-v2"))
    assert d != alternate
    assert d.dataset_id == alternate.dataset_id  # Only the convention ID differs.
    assert d.metadata.basis == alternate.metadata.basis
    np.testing.assert_array_equal(evaluate(d, [0, 1]), evaluate(alternate, [0, 1]))
    assert d.metadata.convention_id == "synthetic-convention-v1"
    assert alternate.metadata.convention_id == "synthetic-convention-v2"
    with pytest.raises(FrozenInstanceError):
        d.metadata.convention_id = "other"


@pytest.mark.parametrize("invalid", ["", "   ", None])
def test_convention_id_must_be_nonempty(invalid):
    with pytest.raises(ValueError, match="convention_id"):
        replace(dataset().metadata, convention_id=invalid)
