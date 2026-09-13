"""Evaluate supplied real nuclear-response data without DM or halo physics.

References guiding the interface (no physical coefficients are bundled):
* Schneck et al., https://arxiv.org/abs/1503.03379v3 (including erratum).
* Fitzpatrick et al., https://arxiv.org/abs/1203.3542v3, Sec. 3 and Appendix A.
* Anand et al., https://arxiv.org/abs/1308.6288, nuclear/WIMP factorization.
* WimPyDD, https://arxiv.org/abs/2106.06207, Appendix A and D.
* Hoferichter et al., https://arxiv.org/abs/1812.05617, Eq. (C14).

Channel strings belong to the declared dataset basis. No implicit translation
between Fitzpatrick responses, proton/neutron form factors, and ChiralEFT4DM
amplitudes is defined. See docs/TARGETS_NUCLEAR_RESPONSES.md before importing.
"""

from dataclasses import dataclass
from types import MappingProxyType
from typing import Literal, Mapping

import numpy as np

from .isotopes import Isotope
from .provenance import SourceReference, _integer, _real, _text


def _array(value, name):
    array = np.asarray(value)
    if array.dtype.kind not in "iuf":
        raise ValueError(f"{name} must be real numerical data (not bool or complex)")
    array = np.asarray(array, dtype=float)
    if np.any(~np.isfinite(array)):
        raise ValueError(f"{name} must be finite")
    return array


def _tuple(value, name):
    array = _array(value, name)
    if array.ndim != 1 or not array.size:
        raise ValueError(f"{name} must be a nonempty one-dimensional sequence")
    return tuple(float(v) for v in array)


def _analytic_domain(value):
    if value is None:
        return (0.0, None)
    endpoints = tuple(value)
    if len(endpoints) != 2:
        raise ValueError("domain requires two endpoints")
    lower = _real(endpoints[0], "domain lower bound")
    upper = None if endpoints[1] is None else _real(endpoints[1], "domain upper bound")
    if lower < 0 or (upper is not None and upper <= lower):
        raise ValueError("domain requires 0 <= lower < upper; use None for an unbounded upper endpoint")
    return lower, upper


def _validated_range(value, domain, source):
    if value is None:
        if source is not None:
            raise ValueError("validated_range_source requires a validated_range")
        return None
    bounds = _tuple(value, "validated_range")
    if (len(bounds) != 2 or not domain[0] <= bounds[0] < bounds[1]
            or (domain[1] is not None and bounds[1] > domain[1])):
        raise ValueError("validated_range must be finite and contained in the mathematical domain")
    if not isinstance(source, SourceReference):
        raise TypeError("validated_range requires its source")
    return bounds


def _in_domain(value, domain, name="domain"):
    value = _array(value, "momentum coordinate")
    if np.any(value < domain[0]) or (domain[1] is not None and np.any(value > domain[1])):
        raise ValueError(f"momentum coordinate outside supplied {name} {domain}")
    return value


def _output(value):
    if np.any(~np.isfinite(value)):
        raise ValueError("response evaluation produced nonfinite values")
    return float(value) if np.ndim(value) == 0 else value


@dataclass(frozen=True)
class MomentumVariable:
    """Map API q² in GeV² to the source coordinate, without using isotope mass.

    ``q_GeV`` means |q|*c; ``q2_GeV2`` means (|q|*c)². In natural units,
    ``y`` = q²*b²/4 (WimPyDD Appendix D) and ``u`` = q²*b²/2
    (Hoferichter et al., Eq. C14). b is always supplied in GeV^-1, with a
    source for its value and any conversion from fm. No empirical b(A) default.
    """

    variable: Literal["q2_GeV2", "q_GeV", "y", "u"]
    b_GeV_inv: float | None = None
    b_source: SourceReference | None = None

    def __post_init__(self):
        if self.variable not in ("q2_GeV2", "q_GeV", "y", "u"):
            raise ValueError("unsupported source momentum variable")
        if self.variable in ("y", "u"):
            b = _real(self.b_GeV_inv, "b_GeV_inv")
            if b <= 0 or not isinstance(self.b_source, SourceReference):
                raise ValueError("y/u require positive b_GeV_inv and its source")
            object.__setattr__(self, "b_GeV_inv", b)
        elif self.b_GeV_inv is not None or self.b_source is not None:
            raise ValueError("b is only applicable to y/u")

    def from_q2(self, q2_GeV2):
        q2 = _array(q2_GeV2, "q2_GeV2")
        if np.any(q2 < 0):
            raise ValueError("q2_GeV2 must be nonnegative")
        with np.errstate(over="ignore", invalid="ignore"):
            if self.variable == "q2_GeV2":
                value = q2.copy()
            elif self.variable == "q_GeV":
                value = np.sqrt(q2)
            else:
                value = (np.sqrt(q2) * (self.b_GeV_inv / 2)) ** 2
                if self.variable == "u":
                    value = 2 * value
        return _output(value)


@dataclass(frozen=True)
class ResponseKey:
    """Ordered channel/isospin key. Integer indices refer to metadata labels.

    For an isoscalar/isovector dataset use labels ('0', '1'); a p/n dataset
    can instead declare ('p', 'n'). No symmetry, signs, or factor two inferred.
    """

    channel: str
    tau: int
    tau_prime: int

    def __post_init__(self):
        _text(self.channel, "channel")
        for name in ("tau", "tau_prime"):
            value = _integer(getattr(self, name), name)
            if value not in (0, 1):
                raise ValueError("isospin indices must be 0 or 1")
            object.__setattr__(self, name, value)


@dataclass(frozen=True)
class ResponseMetadata:
    """Required conventions for raw returned values; none are applied implicitly.

    basis identifies the nuclear-response basis, not a Wilson operator list.
    isospin_definition must spell out operators and coupling normalization.
    normalization must state spin averages, 4*pi, q/m_N powers, and whether
    interference factors are included. source locates those definitions.
    convention_id is a required nonempty, exact-match identifier for this full
    convention (including version). It is not inferred from the descriptive
    fields. Future contractions must check it; no conversion occurs here.
    These are declarations by the importer, not automatically verified physics.
    """

    basis: str
    isospin_labels: tuple[str, str]
    isospin_definition: str
    normalization: str
    response_units: str
    momentum: MomentumVariable
    source: SourceReference
    convention_id: str

    def __post_init__(self):
        for name in ("basis", "isospin_definition", "normalization", "response_units", "convention_id"):
            _text(getattr(self, name), name)
        if isinstance(self.isospin_labels, str):
            raise ValueError("isospin_labels must be two distinct strings")
        labels = tuple(_text(v, "isospin label") for v in self.isospin_labels)
        if len(labels) != 2 or labels[0] == labels[1]:
            raise ValueError("isospin_labels must be two distinct strings")
        if not isinstance(self.momentum, MomentumVariable) or not isinstance(self.source, SourceReference):
            raise TypeError("momentum and source must be immutable metadata records")
        object.__setattr__(self, "isospin_labels", labels)


@dataclass(frozen=True)
class PolynomialResponse:
    """exp(-exponential_decay*x) * sum_n coefficients[n]*x**n.

    x is the dataset's source momentum coordinate. ``domain`` is mathematical:
    None means [0, infinity); (lower, None) has no upper bound. All actual inputs
    must be finite and nonnegative. An optional finite ``validated_range`` is
    separately sourced and always enforced. Its absence is not a validation
    claim. Coefficients use ascending powers. Neither exp(-2*y) nor
    normalization to one is assumed. Signed responses are retained.
    """

    coefficients: tuple[float, ...]
    exponential_decay: float
    domain: tuple[float, float | None] | None
    source: SourceReference
    validated_range: tuple[float, float] | None = None
    validated_range_source: SourceReference | None = None

    def __post_init__(self):
        coefficients = _tuple(self.coefficients, "coefficients")
        domain = _analytic_domain(self.domain)
        decay = _real(self.exponential_decay, "exponential_decay")
        if decay < 0:
            raise ValueError("exponential_decay must be nonnegative")
        if not isinstance(self.source, SourceReference):
            raise TypeError("source must locate polynomial coefficients and domain")
        object.__setattr__(self, "coefficients", coefficients)
        object.__setattr__(self, "domain", domain)
        object.__setattr__(self, "exponential_decay", decay)
        object.__setattr__(self, "validated_range", _validated_range(
            self.validated_range, domain, self.validated_range_source))

    def evaluate(self, x):
        x = _in_domain(x, self.domain)
        if self.validated_range is not None:
            x = _in_domain(x, self.validated_range, "validated_range")
        with np.errstate(over="ignore", invalid="ignore"):
            value = (np.polynomial.polynomial.polyval(x, self.coefficients)
                     * np.exp(-self.exponential_decay * x))
        return _output(value)


@dataclass(frozen=True)
class TabulatedResponse:
    """Piecewise-linear response in the source coordinate; no extrapolation.

    source must locate all nodes and values. Interpolation is a numerical
    approximation; its physical accuracy is the data provider's responsibility.
    The finite node domain is always enforced. If supplied, a separately sourced
    finite validated_range inside that domain is also always enforced.
    """

    nodes: tuple[float, ...]
    values: tuple[float, ...]
    source: SourceReference
    validated_range: tuple[float, float] | None = None
    validated_range_source: SourceReference | None = None

    def __post_init__(self):
        nodes, values = _tuple(self.nodes, "nodes"), _tuple(self.values, "values")
        if len(nodes) < 2 or len(nodes) != len(values):
            raise ValueError("at least two matching nodes and values required")
        if nodes[0] < 0 or any(b <= a for a, b in zip(nodes, nodes[1:])):
            raise ValueError("nodes must be nonnegative and strictly increasing")
        if not isinstance(self.source, SourceReference):
            raise TypeError("source must locate tabulated nodes and values")
        object.__setattr__(self, "nodes", nodes)
        object.__setattr__(self, "values", values)
        object.__setattr__(self, "validated_range", _validated_range(
            self.validated_range, self.domain, self.validated_range_source))

    @property
    def domain(self):
        return self.nodes[0], self.nodes[-1]

    def evaluate(self, x):
        x = _in_domain(x, self.domain)
        if self.validated_range is not None:
            x = _in_domain(x, self.validated_range, "validated_range")
        return _output(np.interp(x.ravel(), self.nodes, self.values).reshape(x.shape))


@dataclass(frozen=True)
class NuclearResponseDataset:
    """One isotope and convention, arbitrarily many explicitly supplied entries.

    Multiple datasets for the same isotope coexist as independent objects.
    A missing entry raises KeyError (including an absent transposed entry).
    Supply an explicitly sourced zero polynomial for a known vanishing entry.
    There is no isotope-composition weighting and no global dataset registry.
    """

    dataset_id: str
    isotope: Isotope
    metadata: ResponseMetadata
    responses: Mapping[ResponseKey, PolynomialResponse | TabulatedResponse]

    def __post_init__(self):
        _text(self.dataset_id, "dataset_id")
        if not isinstance(self.isotope, Isotope) or not isinstance(self.metadata, ResponseMetadata):
            raise TypeError("isotope and metadata must be immutable records")
        responses = dict(self.responses)
        for key, response in responses.items():
            if (not isinstance(key, ResponseKey)
                    or not isinstance(response, (PolynomialResponse, TabulatedResponse))):
                raise TypeError("responses must map ResponseKey to immutable polynomial/table data")
        object.__setattr__(self, "responses", MappingProxyType(responses))

    def evaluate(self, *, channel, tau, tau_prime, q2_GeV2):
        """Return raw W for this dataset's isotope, preserving the q² input shape."""
        key = ResponseKey(channel, tau, tau_prime)
        try:
            response = self.responses[key]
        except KeyError:
            raise KeyError(f"{self.dataset_id}: missing response {key}") from None
        return response.evaluate(self.metadata.momentum.from_q2(q2_GeV2))
