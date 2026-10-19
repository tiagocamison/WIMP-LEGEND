# CODE-AUDIT-003 detector remediation

This is a coder remediation record for `CODE_AUDIT_003_GENERIC_ANALYSIS.md`,
not an independent acceptance verdict. The repair branch is based directly on
audited commit `166b1855ce75600bafb7a39947bd0887c31d8617`. The audit's two MAJOR
findings were reproduced: seven unchanged external cases failed before repair.
The public API is unchanged. No physical isotope/response or detector data is
introduced, and no rate, interaction, statistics or sensitivity code is changed.

## GA-AUDIT-001: propagate observed efficiency boundaries

`DetectorResponse.expected_counts` now combines declared observed-efficiency
boundaries inside the observed window with its bin edges before mapping them
back through `ConstantQuenching` or `LindhardQuenching`. Deterministic responses
split the truth integral at those crossings. Constant Gaussian widths add
crossings at the same boundaries and their +/-1,2,4,8 sigma neighborhoods.
These are integration hints, not a truncation or normalization of the Gaussian.

The unchanged audit cases recover **0.01 events** for flat unit truth on [0,1]
with acceptance only on (0.4,0.41), using NoSmearing, sigma=0, and sigma=1e-5.
Callers of custom quenchers, callable widths or custom kernels still need to
provide the truth breakpoints relevant to those models explicitly.

## GA-AUDIT-002: use local Gaussian probability coordinates

`GaussianResolution.bin_probabilities` splits each efficiency-weighted bin at
declared efficiency boundaries and the Gaussian mean in observed energy,
before transforming any interval. Each segment uses CDF probabilities on the
negative side of the mean or survival probabilities on the positive side.
The inverse therefore receives probabilities at most one half instead of a
coarse-bin CDF that rounds to one. Segment contributions are summed with
`math.fsum`. The unweighted and ConstantEfficiency shortcuts remain unchanged.

If local probability interpolation underflows while the segment mass is still
representable, the inverse uses log probabilities, including the lower endpoint
log CDF. Nonfinite reconstructed coordinates raise explicitly. Neither
probabilities nor callback results are clipped. Existing finite/probability
validation and IntegrationWarning-as-error behavior are preserved.

The unchanged coarse-bin audit cases retain accepted probabilities
**7.619661958202802e-24**, **6.219831985865741e-16**, and
**1.2791904478283997e-12** for means 0.30,0.32,0.33 with sigma=0.01.
The full fold with explicit truth hints also recovers **0.01 events**.
This remains a floating-point calculation: extremely remote normal tails can
underflow to zero, as previously disclosed. Undeclared narrow callback features
are not guaranteed to be discovered by adaptive integration.

## Regression evidence

`tests/test_detector_audit_regressions.py` adds 20 cases. Its references use
independent `math.erf`/`erfc` normal integrals, finite Poisson sums, and explicit
likelihood ratios. No expected value calls the production detector or statistics
function under test. Coverage includes the original seven counterexamples,
left/right tail windows, finite in-bin callback coordinates, compatible bin
partitioning, duplicate/outside hints and caller ownership, both built-in
quenchers, and a nonunit-quenching narrow Gaussian.

The discontinuous-efficiency pipeline has truth density 10 events/(kg day keV),
truth range [0,1], exposure 10 kg day, unit quenching, sigma=0.01 keV,
acceptance on (0.4,0.41), and bins [0,0.405,1]. Its independently known counts
are **S=(0.5,0.5)** and **B=(0.81,1.19)** with flat background density
0.2 events/(kg day keV_obs); observed counts are (1,1). Observed and expected
limits from both backends are compared with separate mathematical references.
The entire pipeline is repeated with epsrel=1e-10 and root xtol/rtol=1e-13,
checking both intermediate counts and final limit convergence.

| Independent limit reference | mu upper |
| --- | ---: |
| Exact counting, observed and discrete median expected | 3.877404490410548 |
| Asymptotic shape, observed | 2.4479166828589265 |
| Asymptotic shape, Asimov expected | 2.3766484464804964 |

All commands below ran with warnings treated as errors. The test dependencies
were supplied through a temporary PYTHONPATH; package dependencies were not
changed. Environment: Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0, pytest 9.1.1,
and mpmath 1.4.1 for the external oracle only.

| Validation | Result |
| --- | --- |
| Exact seven external counterexamples before repair | 7 failed, 117 deselected |
| Exact same counterexamples after initial repair | 7 passed, 117 deselected |
| Targeted detector/new-regression/sensitivity suite | 97 passed |
| `python -m pytest -q -W error` (entire repository) | **675 passed** (655 inherited + 20 added) |
| Unchanged pre-rate and rate external suites | **391 passed** (154 + 237) |
| Unchanged generic audit Appendix A, including the smooth full pipeline and all seven counterexamples | **124 passed** |

The full repository and external suites therefore executed **1,190 passing
cases** in total. The targeted run is a subset and is not added to that total.
No inherited repository test or external assertion/tolerance was altered.
`git diff --check` also passes. The frozen NREFT normative document remains
byte-identical (SHA-256
`a2c57658add45b273ca082d8b01fbab7dd739f8625d6f04a438d619d4cf385cd`).

## External test identity

All inherited tests and tolerances are retained. The audit's Appendix A is
extracted verbatim outside the production repository, including its high-
precision mpmath oracle (80 decimal digits). No audit assertion is skipped or
marked xfail. SHA-256 identities:

| External script | SHA-256 |
| --- | --- |
| Pre-rate independent | `1f1827d261b49a0928db7eec5c33c0bc5d518d976d107f245fd4de6639ac23e6` |
| Pre-rate falsifications | `99438a0db818631638f99d80b7282c8f6ed8d56505443abf4057b4e0005de6cb` |
| Rate audit | `a90623aca9133d9293181cad00c4db5dbbb5fbc0f72cc72793612c172ba24aae` |
| Generic audit Appendix A | `1705ebe6e24403a54573d3dd6b0ba8724c681bbc29a1e684a93e8bdc1fc05aec` |

## Statistical conditions remain explicit

GA-AUDIT-003 is retained as a scope condition: the shape backend is asymptotic,
has `status='asymptotic_unvalidated'`, and does not certify sparse-count coverage.
The added pipeline explicitly checks that metadata. No new low-count shape
calibration is claimed. The audit's coverage study is not rerun because the
statistics implementation and calibration are unchanged.

GA-AUDIT-004 is also retained: shape-root absolute tolerance is in mu, whereas
counting solves in signal events. Template normalization and explicit root
tolerances require convergence checks for tiny mu. The unchanged external
template-rescaling test with tighter root tolerances is included in validation.

An independent re-audit is still required for milestone acceptance. A green
generic suite does not validate nuclear data, experimental calibration,
composition, backgrounds, uncertainty models or a physical Ge analysis.
