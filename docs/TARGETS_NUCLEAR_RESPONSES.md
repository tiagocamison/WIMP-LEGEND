# Targets and external nuclear responses

This layer contains **no production isotope database or physical response coefficients**.
It validates and evaluates caller-supplied data, preserving their provenance. All
numerical fixtures in the new tests are synthetic, including their isotope labels,
masses, spins, fractions, oscillator lengths, coefficients, and domains.

## Public API and data flow

Exports from `wimp_legend.targets`:

| API | Responsibility |
| --- | --- |
| `SourceReference(reference, locator, kind, upstream=(), notes="")` | Immediate source, precise location, upstream calculation, unresolved qualifications |
| `Isotope(symbol, A, Z, mass_GeV, spin, sources)` | One ground-state nuclide; separate source entries for `identity`, `mass_GeV`, `spin` |
| `TargetComposition(components, fraction_kind, source)` | Ordered `(isotope, fraction)` pairs; explicit `number` or `mass` |
| `TargetComposition.from_reported(components, fraction_kind, source, *, normalize=False)` | Strict import by default; explicit normalization with `normalize=True` |
| `composition.normalized()` | New immutable normalized composition with derivation provenance |
| `to_mass_fractions()`, `to_number_fractions()` | Explicit conversions, returning immutable compositions |
| `MomentumVariable(variable, b_GeV_inv=None, b_source=None)` | Source coordinate and any oscillator length with provenance |
| `ResponseKey(channel, tau, tau_prime)` | Ordered, dataset-local response entry |
| `ResponseMetadata(basis, isospin_labels, isospin_definition, normalization, response_units, momentum, source, convention_id)` | Required convention declarations and exact-match machine identifier |
| `PolynomialResponse(coefficients, exponential_decay, domain, source, validated_range=None, validated_range_source=None)` | Analytic response with separate mathematical domain and source-supported range |
| `TabulatedResponse(nodes, values, source, validated_range=None, validated_range_source=None)` | Finite interpolation domain, optionally restricted by a source-supported range |
| `NuclearResponseDataset(dataset_id, isotope, metadata, responses)` | One isotope and convention; independently selectable calculation |
| `dataset.evaluate(channel=..., tau=..., tau_prime=..., q2_GeV2=...)` | Raw response, scalar or shape-preserving NumPy array |

External property data feed `Isotope`; external sample fractions feed composition.
External nuclear-structure calculations feed datasets. A future rate calculation
will select a dataset for each isotope, combine its responses with separately
computed DM responses and halo moments, then sum isotope rates using the appropriate
fractions. None of those contractions or rates are implemented here.

Targets import only local provenance helpers, the standard library, and NumPy.
They do not import halo, kinematics, interactions, rates, or detector code. The
linear dependency sketch in the older `ARCHITECTURE.md` should not be interpreted
as an import chain: these physics inputs meet at the future rate calculation.

## Frozen interface conventions

- `mass_GeV` means bare-nucleus rest energy in GeV; `spin` means nuclear J in
  units of hbar. No atomic mass conversion or mass approximation is automatic.
- Ground-state identity is `(Z, A)`. Duplicate detection uses this identity even
  if masses, sources, symbols, or spins differ. Isomers are not represented.
- Validation is structural, not a periodic-table or nuclear-database lookup.
  Importers must verify symbol/Z consistency against their external source;
  compositions also reject contradictory symbol/Z assignments internally.
- Fractions describe the nuclear ensemble. For number fractions x and nuclear
  mass fractions w: `w_i=x_i*m_i/sum(x*m)` and
  `x_i=(w_i/m_i)/sum(w/m)`. Conversion uses the supplied masses, never A.
  Equating these nuclear mass fractions with measured bulk detector-mass fractions
  requires a later explicit electron/binding convention.
- Fractions must be finite, between zero and one, and sum to one within absolute
  `1e-12` with zero relative tolerance. This is a numerical policy, not a physical
  constant. Accepted input values remain unchanged; no implicit renormalization.
  Explicit conversions normalize computed weights and create derived source
  records linked to the original. See the composition audit trail below.
- Input is always `q2_GeV2=(|q|c)^2`, not recoil energy or signed momentum.
  Scalar input returns Python float; arrays preserve arbitrary shape.
- The source coordinate is `q2_GeV2`, `q_GeV`, `y=q²b²/4`, or `u=q²b²/2`.
  The latter definitions follow [WimPyDD, Appendix D](https://arxiv.org/pdf/2106.06207)
  and [Hoferichter et al., Eqs. C14–C15](https://arxiv.org/pdf/1812.05617),
  respectively. b is supplied in GeV^-1 with a source; there is no b(A) default.
- Polynomial evaluation is `exp(-decay*x)*sum(a_n*x**n)`; coefficients are
  ascending powers. `decay` is mandatory, with reciprocal-x units when x has
  units. Source locators must cover every coefficient and decay. A mathematical
  domain must not be misrepresented as a source-validated range.
  Tables interpolate linearly in x, not necessarily in q².
- Analytic `domain=None` means mathematical `[0, infinity)`; `(lower, None)`
  permits a nonzero lower bound without an upper bound. Explicit finite domains
  remain supported. Use `None`, not floating-point infinity, for unboundedness.
  Evaluated inputs must always be finite and nonnegative. Absence of a cutoff is
  not evidence of physical validity at arbitrarily large momentum.
- Tables retain their finite interpolation domain `[nodes[0], nodes[-1]]` and
  forbid extrapolation. Both response types optionally accept a finite
  `validated_range` contained in their mathematical/interpolation domain, in
  the same source-coordinate units. It requires `validated_range_source` and
  is **always enforced**, not advisory or configurable. Its absence means no
  validated range has been recorded; the code does not invent one. Source
  records are declarations, not independent certification.
- Out-of-domain or out-of-validated-range evaluation raises `ValueError` rather
  than clamping or returning zero. Nonfinite results are rejected; extreme
  polynomial/exponential cancellation is not specially resummed.
- Responses are real and may be signed. There is no clipping, spin averaging,
  normalization at zero momentum, abundance weighting, or automatic q/m_N factor.
- Indices 0 and 1 select the two **declared labels**. Labels `('0','1')` can
  denote isoscalar/isovector; `('p','n')` denotes a different convention. The full
  definition and coupling normalization must be supplied. No basis conversion.
- Keys are ordered. Neither `(0,1)` versus `(1,0)` nor reversed interference
  channel names are identified. Missing data raise `KeyError`; a known zero
  needs an explicitly sourced zero polynomial. Empty datasets are permitted.
- Dataset basis identifiers are explicit strings, with no implicit aliases.
  The additional required `convention_id` is immutable, nonempty and compared
  exactly. It identifies the entire versioned normalization/isospin convention,
  not merely a channel vocabulary. Evaluation preserves it without conversion.
  A future R/W contraction must require an exact convention-ID match; a shared
  basis name or channel name is insufficient. No canonical physical ID or
  matching adapter is introduced in this pass.
  Suggested ASCII names for the conventional six channels are `M`,
  `Sigma_prime`, `Sigma_double_prime`, `Phi_double_prime`, `Phi_tilde_prime`,
  `Delta`; ordered interference examples are `Phi_double_prime_M` and
  `Delta_Sigma_prime`. Availability and definitions remain dataset-specific.
- Lists/arrays become owned tuples; mappings become read-only copies containing
  immutable records. Arbitrary callables are not accepted as response data.
  Frozen objects are not all hashable; use isotope identity/dataset_id for lookup.

## Explicit composition normalization and derivations

The strict constructor still rejects rounded totals such as 0.9999 or 1.0001.
Use `TargetComposition.from_reported(pairs, kind, source, normalize=True)` to
normalize those reported inputs. Without that explicit flag the factory is
strict too. All component validation runs before normalization: no negatives,
nonfinite fractions, fractions above one, duplicate isotopes, invalid fraction
kind, empty input, or zero total. The explicit operation accepts any positive
total; it does not impose an undocumented rounding tolerance. Inputs remain
fractions, not percentages or arbitrary weights above one.

For an existing valid composition, `normalized()` returns a new object with the
same fraction kind. It records the operation even when the original sum is
exactly one. Constructor-accepted roundoff can therefore be normalized explicitly
without mutating the original object. Conversion to the already-current fraction
kind still returns the existing object without recording a fictitious conversion.

`SourceReference` gains three optional immutable fields; existing construction
remains compatible:

| Field | Meaning |
| --- | --- |
| `derived_from` | The complete preceding SourceReference, including its original citations/notes |
| `transformation` | `None` for directly supplied values, or `normalize_number`, `normalize_mass`, `number_to_mass`, `mass_to_number` |
| `input_fractions` | Owned tuple of the values **before** this operation, in unchanged isotope order |

Derived records have kind `implementation` and locate the internal transformation,
rather than attributing computed numbers to the original publication. The previous
record is retained once through `derived_from`; its bibliographic `upstream`
citations are not flattened or duplicated. Each operation retains the isotope
objects, including their masses and property sources. A chain of normalization,
number-to-mass conversion and mass-to-number conversion thus preserves the original
rounded inputs and every intermediate input, with unambiguous operation tags.
No general-purpose graph machinery or new top-level model class is added.

For example, with caller-supplied isotope objects and a source record:

```python
reported_pairs = [(isotope_a, 0.25), (isotope_b, 0.7499)]  # synthetic example
number = TargetComposition.from_reported(
    reported_pairs, "number", reported_source, normalize=True,
)
mass = number.to_mass_fractions()
assert mass.source.transformation == "number_to_mass"
assert mass.source.derived_from is number.source
assert number.source.input_fractions == (0.25, 0.7499)
assert number.source.derived_from is reported_source
```

## Reference audit and consequences

References were consulted on 2026-09-12. “First located” below describes this
audit's source trail, not a claim of historical priority. No external library code
or physical coefficient arrays were copied into executable package code.

1. [Schneck et al., arXiv:1503.03379v3](https://arxiv.org/abs/1503.03379v3):
   latest version, including the erratum. Used for operator/channel terminology,
   isospin definitions, and interference. The erratum reports a MATLAB form-factor
   bug affecting interference, checked against Mathematica. Section I uses
   `c^0=(c^p+c^n)/2`, `c^1=(c^p-c^n)/2`. Its references 21–23 lead to the
   foundational response calculations; it is not our coefficient database.
2. [WimPyDD paper](https://arxiv.org/pdf/2106.06207) replaces the Hepforge page,
   which returned a bot challenge. Sections 2–3 and Appendix A separate Wilson,
   nuclear, and halo inputs; cached experimental responses are distinct from W.
   Appendix D gives isotope-specific polynomial inputs. Its section 3.5.2
   example uses `c_tau=cp*[1+r,1-r]`, differing from Schneck's half-sum definition.
   This motivates explicit convention metadata, not an automatic conversion.
   References 10 and 23 identify Anand et al. and Catena–Schwabe as response sources.
3. [ChiralEFT4DM project page](https://theorie.ikp.physik.tu-darmstadt.de/strongint/ChiralEFT4DM.html)
   identifies coherent one-/two-nucleon chiral-EFT responses and quark/gluon
   matching, citing PRD 99, 055031; PRL 122, 071301; PRD 94, 063505; PLB 746, 410.
   [Hoferichter et al., arXiv:1812.05617v2](https://arxiv.org/abs/1812.05617v2)
   provides the primary calculation and the `F_M^±`, `F_Phi''^±`, `F_pi`, `F_b`
   structure factors. These are not interchangeable with generic bilinear W
   entries. Merely changing the momentum coordinate does not convert an amplitude
   into a response. A future chiral adapter needs its own basis and matching rules;
   two-nucleon channels must not be forced into fictitious tau pairs here.
4. [Pinned RAPIDD tree](https://github.com/cheekyparticle/RAPIDD_for_DM/tree/02035a29e4535e42ec21363f369f5ea431bc4482)
   is a historical implementation/surrogate framework. Its
   [effFormFact.c](https://github.com/cheekyparticle/RAPIDD_for_DM/blob/02035a29e4535e42ec21363f369f5ea431bc4482/lib/source/effFormFact.c)
   header cites 1203.3542, 1308.6288 and 1605.08043. A file-wide citation does not
   prove every branch shares a normalization or establish every coefficient's origin.

The foundational [Anand et al. preprint](https://arxiv.org/pdf/1308.6288) also
needs care: Eqs. 69–70 and the nearby `SetCoeffsNonrel` example appear inconsistent
about a factor of two. The published PRC 89, 065501 and the actual versioned
Mathematica code must be reconciled before implementing a convention adapter.
Its b(A) default appears in the `SetIsotope` documentation, but is not an exact
isotope property. [Fitzpatrick v3](https://arxiv.org/abs/1203.3542v3) explicitly
notes corrections to factors of 4*pi, another reason to pin versions.

## Notebook quantities and coefficient trail

Notebook cell indices below are zero-based. Source locators cover each member
of the listed coefficient families; numbers in this section are an audit record,
not importable or trusted package data.

| Quantity first located in notebook | Classification and convention | Trace and disposition |
| --- | --- | --- |
| Cell 6: `A=76`, `76*0.9315`; cell 12: `Z=32` | Identity plus approximate mass model, GeV | Identity is Ge-76; mass is an approximation, not a sourced nuclear mass. No production value supplied. |
| Cell 25: `75.9214*0.9315` | Unsourced isotope-mass input | Looks like an atomic-mass-based value, but no mass dataset, uncertainty, or atomic/nuclear conversion is specified. Unresolved. |
| Cell 12: `Z*0.938+(A-Z)*0.940` | Additive free-nucleon mass approximation in GeV | RAPIDD `f`, line 132, uses this formula; it is not the bound-nucleus mass. No exact-mass authority traced. |
| Cell 12: `41.467`, `45`, `25`, `0.1973` in b(A) | Empirical oscillator prescription and rounded unit conversion | Anand `SetIsotope` documentation locates the b(A) formula; RAPIDD `f` has a commented `0.197` conversion and active B2 lookup. Notebook uses `0.1973`. These are not interchangeable precision conventions. No default encoded. |
| Cell 12: pp, pn, nn arrays below | Published-response-like p/n input, `exp(-2*y)` | Exact notebook arrays found in RAPIDD `FF76`, lines 898–904. Primary trail reaches Fitzpatrick Appendix A.3, Eq. 81, with constant-term differences described below. Not admitted. |
| Cell 12: `pp+2*pn+nn`, divided by `76**2=5776` | Equal-p/n SI coupling and form-factor normalization assumption | Notebook calculation, not raw W. The factor two and A² normalization belong to a later documented contraction. |
| Cell 21: `cp`, `cn` and their selected entries | Particle-model assumptions | Future Wilson input / NREFT engine; not isotope data or nuclear-structure coefficients. |
| Cells labeled `enrGe` / `natGe` | Background configurations | They do not supply a complete isotope composition. Natural/enriched fractions remain user inputs. |

Each row below lists **all** coefficients in ascending powers, including appended
zeros, as observed in notebook cell 12 and pinned RAPIDD `FF76`:

| Entry | Coefficients a0 through a10 |
| --- | --- |
| pp | 1024, -2800, 2900, -1400, 340, -40, 1.8, -0.0024, 0, 0, 0 |
| pn | 1408, -4200, 4800, -2600, 730, -100, 6.5, -0.11, 0, 0, 0 |
| nn | 1936, -6300, 7800, -4700, 1500, -260, 21, -0.66, 0.0069, 0, 0 |

[Fitzpatrick Eq. 81](https://arxiv.org/pdf/1203.3542v3) gives the same
nonconstant terms, but rounded constants 1000, 1400, 1900. RAPIDD/notebook constants
equal Z², ZN, N² for Z=32, N=44. Interpreting that as restoration of exact
zero-momentum counting is an inference, not a documented derivation in this audit.
Section 3.3 describes the underlying shell-model calculation; Eq. 81 uses F in a
p/n basis, not directly this API's desired isoscalar/isovector W convention.
Coefficient provenance is therefore partly traced, normalization/import remains open.

Other inspected numerical families are catalogued without promotion:

| Source/location | Encountered quantities | Authority/convention/original calculation |
| --- | --- | --- |
| RAPIDD `effFormFact.c`, B/B2 declarations at lines 47/49 | Oscillator lookup arrays, indexed by A-1 | Implementation input; related commented b(A) prescription identified, exact table-generation provenance not established. |
| Same file, `FF70`, `FF72`, `FF73`, `FF74` M-pp entries | Neighboring Ge polynomial coefficient lists | Implementation in p/n convention; header-level primary citations only checked for these entries. No individual validation or import. |
| Same file, `FF76`, following M block | Ppp and M/Ppp coefficient families | Historical Phi'' / interference labels; potential primary Eq. 81 trail, but ordering/normalization not independently reconciled. Excluded. |
| WimPyDD Appendix D, Xe.tab example | Xe isotope A, fractional abundances, spins, internal itar identifiers | Implementation example, not an authoritative abundance/mass evaluation. Original property dataset not identified there. Not imported. |
| Schneck discussion of nuclear expectation values | Illustrative spin/orbit response magnitudes | Secondary discussion pointing to Fitzpatrick 1211.2818; no numerical input extracted for this layer. |
| Hoferichter Table VI | Maximal polynomial powers by element and channel | Primary fit specification associated with C14; no numerical fit imported. |

The wider paper tables and unrelated repository rate/background constants have
not been audited as a nuclear-property database. No claim of exhaustive external
coefficient validation is made. In particular, WIMpy in the notebook and WimPyDD
in the supplied reference are distinct projects; agreement with one does not
validate the other's conventions.

## Required decisions before the NREFT engine

1. Select versioned nuclear-mass and spin data, including atomic-to-nuclear
   conversion constants and uncertainties. Supply the full sample composition
   with original number/mass semantics; do not invent LEGEND enrichment.
2. Choose the canonical Hamiltonian/isospin convention, spelling out t0, t1,
   proton sign and the full/half-sum mapping. Resolve the literature/code ambiguity.
3. Choose raw W normalization: reduced matrix elements, spin averages, 4*pi,
   momentum factors outside versus inside W/R, and signed ordered interference
   definitions. Verify it against a primary calculation and one independent code.
4. Pin each dataset's files/version, source equations, isotope, oscillator length,
   coefficient precision, any source-supported validated range, and
   missing-versus-exact-zero entries. Leave unpublished cutoffs unspecified.
   Validate a physical dataset separately before using it for limits.
5. Freeze the DM operator list (including O2 policy), DM spin, momentum sign,
   nucleon reference mass, Wilson units, and any matching scale. None are target
   assumptions and none are implemented here.
6. Define composition-to-rate factors per kg and recoil-unit conversions once;
   avoid double-counting isotope abundances or nuclear spin factors.
7. Keep chiral amplitudes/two-body responses and generic one-body W datasets
   distinct. Design and validate explicit matching adapters when required.

## Tests

`tests/test_isotopes.py` checks property validation, per-property provenance,
normalization, duplicates, fraction semantics, analytic/randomized conversion
round trips, multiple elements, zero fractions, and immutable input ownership.
`tests/test_nuclear_responses.py` checks polynomial and table evaluation, all
coordinate definitions, an independent q²/y reparametrization, array shapes,
signed responses, missing/ordered entries, invalid indices, strict domains,
metadata preservation, mutation resistance, and composition independence.
No test asserts that the notebook's physical coefficients are correct.

Hardening validation: **191 tests passed** with warnings treated as errors.
Command: `python -m pytest -q -W error` (runtime test dependencies supplied via
PYTHONPATH). This adds 33 cases to the previously passing 158-test suite.

The only alterations to pre-existing tests are:

- The shared response-dataset fixture now supplies the newly required synthetic
  `convention_id`; all tests using it retain their original assertions.
- The fraction-conversion test now checks `m.source.derived_from is c.source`
  and the `number_to_mass` tag, replacing the obsolete expectation that converted
  values reuse the original source object as their immediate source.

No old tests were removed. The invalid-infinite-domain-endpoint test is retained:
unbounded analytic endpoints use `None`, while explicit numeric endpoints stay
finite. New tests cover unbounded/negative coordinates, enforced and separately
sourced validated ranges, finite tables, convention-ID identity/immutability,
rounded normalization with both fraction kinds, unchanged inputs, exact sums,
invalid normalization inputs, and the full derivation chain in both directions.

Before an O1 contraction, the remaining blocker is a physically verified, versioned
normalization/isospin contract and its canonical convention ID, paired with a
traceable nuclear dataset and an independent normalization check. The prior
Anand/Fitzpatrick/WimPyDD ambiguity remains unresolved. No physical coefficients,
R functions, Wilson coefficients, or automatic convention conversions were added.
