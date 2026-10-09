# GENERIC-ANALYSIS-001 implementation report

**Coder implementation and verification; not independent certification.**
Date: 2026-10-09. Repository: https://github.com/tiagocamison/WIMP-LEGEND.

## Source identity and staged history

The restored checkout's `refactor` and `origin/refactor` tips both exactly match
**44c4f2523b83b926a8b11d153ebfd5d8d04531a5**, parent
**efd56745fa5e35461122ae36afcfe105115061b1**. This is the exact audited base,
not a newer unpinned substitute. Initial tracked/untracked repository status was
clean. Work is isolated on **generic-analysis-001**, directly based on that SHA.
No main/refactor mutation, push, merge or history rewrite is performed.

Required architecture/rate/convention/resolution documents, relevant tests and
all existing analysis skeletons were inspected before editing. The two audit
reports were recovered from the supplied file collection. Their original test
blocks were extracted outside the checkout and match the published SHA-256
hashes. The existing detector/background/statistics tests and skeletons were
empty. No nonempty inherited test was modified.

Separate commits:

1. **be9825f7350865bcb586934c50bd2c2e11b27532** — Phase-A rate range fix,
   permanent reproducer, green regression gate and resolution record.
2. **10fcb541575ba43c3a8ce01a7eef8087c9a37c56** — detector ports/models and
   generic numerical/count validation utility, with 54 initial detector cases.
3. **e6fc884a2cd9c89f36f1726c1a332d83d57507d2** — backgrounds and statistical
   backends, with 15 background and 59 statistics cases.
4. Final commit — sensitivity integration, scan metadata, narrow-window Gaussian
   integration hints, eight additional detector cases, end-to-end tests, example
   and documentation. Its exact SHA is supplied in the delivery message; this
   report is itself included in that commit. `git rev-parse HEAD` identifies it
   after checkout of the delivered branch/patch series.

Unlike the squashed history discussed in RATE-AUDIT-002, this branch preserves
an actual separately committed green Phase-A checkpoint before Phase B.

## Phase A: original defect and closure

Before edits: **492 repository tests pass**; original independent audit scripts
**154 pass**, second rate audit **236 pass / 1 fails**. The failing original
RATE-AUDIT-001 reproducer returns zero instead of
**2.669722691466525e-291 events/(kg day keV)**. The area conversion erased a tiny
natural-unit integral before the larger target-count factor could restore it.

The final product uses frexp/ldexp mantissa/exponent arithmetic, with exactly
these factors: natural speed-flux integral, (hbar c)^2, 1e-6 recoil Jacobian,
c in cm/s, rho/m_chi, nuclei/kg, seconds/day. The mass division is incorporated
separately so the density/mass ratio does not itself overflow/underflow first.
Zero and signs are preserved. Genuine overflow and a nonzero result rounding to
zero raise ValueError. Representable subnormals retain ordinary reduced relative
precision. This does not fix upstream response/callback/quadrature underflow,
ill-conditioned cancellation, or arbitrary-precision errors.

Verified corrected output: **2.6697226914665258e-291**, within the original
90-digit Decimal oracle's 1e-12 relative tolerance (absolute tolerance zero).
The permanent original reproducer is copied verbatim, with its synthetic fixture,
into tests/test_rate_product.py; additional scaling/sign/Decimal range tests are
included. The original external scripts remain unchanged.

At the Phase-A gate: **504 repository tests pass**, **154 old audit tests pass**,
**237 rate-audit tests pass**, all `-W error`, before commit be9825f and before
any Phase-B code. See RATE_AUDIT_002_RESOLUTION.md for commands and hashes.

## Public API inventory and data flow

Full signatures, units, semantics and examples are in the linked documents:

| Layer | Public exports / operations | Detailed contract |
| --- | --- | --- |
| Detector | IntegrationPolicy; QuenchingModel, ConstantQuenching, LindhardQuenching; ResolutionModel, NoSmearing, GaussianResolution; Efficiency, ConstantEfficiency, CallableEfficiency; CountResponse, DetectorResponse, KernelResponse, MatrixResponse | [DETECTOR_RESPONSE.md](DETECTOR_RESPONSE.md) |
| Detector methods | mean_energy_keV, quenching_factor, sigma_at, bin_probabilities, expected_counts, apply_truth_counts | truth/mean/reconstructed energies explicit; expected-count vector output |
| Background | BackgroundModel, ZeroBackground, FlatBackground, CallableBackground, BinnedBackground, CompositeBackground; expected_counts, differential_rate | [BACKGROUNDS.md](BACKGROUNDS.md) |
| Statistics | StatisticsBackend, RootPolicy, LimitResult, FitResult, PoissonCLsCounting, BinnedPoissonProfile | [STATISTICS.md](STATISTICS.md) |
| Statistics methods | observed_upper_limit, expected_upper_limit, acceptance_margin, acceptance_metadata; log_cls; log_likelihood, fit, q_mu | exact counting and asymptotic shape independently named |
| Sensitivity | SensitivityAnalysis, AnalysisResult, GridResult, ParameterScan, scan_parameter | [SENSITIVITY.md](SENSITIVITY.md) |
| Sensitivity methods | observed_limit, expected_limit, expected_grid; accepted_sample_ranges | delegates models and count-level inference |
| Existing rate | differential_rate_per_kg_day_keV API unchanged | final product arithmetic only |

Dependency direction:

- Detector/background/statistics depend on generic numerical/count utilities,
  their own layer and NumPy/SciPy only.
- Sensitivity depends on those ports and supplied callables/counts, without
  importing rate/halo/target/interaction physics.
- The rate/particle/foundation layers never import the new analysis layers.
- Custom pure models/backends need no subclassing or mutable global registry.

A shared output is an owned 1D ndarray of **expected events per observed bin**;
immutable result records retain tuples. No general-purpose hierarchy, empirical
registry, hidden mutable state or new numerical dependency is introduced.

## Equations, units and frozen choices

Input truth spectrum r: **events/(kg day keV_nr)**. Exposure X: **kg day**.
E_nr is true nuclear keV, E_ee is quenched mean electron-equivalent keV, E_obs is
reconstructed keV. Truth range is explicit and nonnegative; finite negative
observed edges are legitimate Gaussian outcomes.

\[
S_j=X\int r(E_{\rm nr})\epsilon_{\rm truth}(E_{\rm nr})
\int_{a_j}^{b_j}K(E_{\rm obs}\mid E_{\rm nr})\epsilon_{\rm obs}(E_{\rm obs})
\,dE_{\rm obs}\,dE_{\rm nr}.
\]

Constant q is explicit in (0,1]; zero quenching is unsupported. Lindhard Z and k
are explicit; the formula specified in the task is implemented with zero-recoil
handling and no Ge k. Gaussian sigma is a constant or function of **E_ee_mean**,
not E_nr or E_obs. Bin probabilities use CDF/survival differences, never window
renormalization. Arbitrary observed efficiency is integrated inside each bin.
Efficiencies are strictly [0,1]. Deterministic bins are [a,b), final upper edge
included. Known monotone quenchers split bin crossings; constant Gaussian widths
also supply scale-aware truth breakpoints, without truncating their tails.

Custom kernel/matrix efficiency ownership is required explicitly. An embedded
acceptance rejects an additional truth efficiency. Separate observed efficiency
requires sub-bin information and belongs inside the supplied kernel/matrix.
Matrix M[j,i] gives observed=M@truth_counts, column probabilities <=1 (8 eps
sum roundoff tolerance, no renormalization). Its within-bin averaging description
is required, and its exact edges/full truth range are enforced.

Background differential rate b: **events/(kg day keV_obs)**;
B_j=X integral_bin b dE_obs. Binned mode is explicitly per_unit_exposure or
fixed_expected_counts at a required reference exposure. Fixed counts do not
silently scale, and incompatible edges do not silently rebin.

Statistics uses known finite nonnegative S/B, actual integer n, and
lambda_i=mu*S_i+B_i. Confidence default .90 is recorded explicitly.
Counting CLs=P(N<=n|B+mu*S)/P(N<=n|B), with root CLs=1-CL. It collapses bins,
uses stable log-CDF/series arithmetic, and solves in signal events. Expected
median uses the smallest integer Poisson .5 quantile, not n=B.

Shape log L=sum(n log lambda-lambda-log(n!)); fit mu_hat>=0;
q_mu=0 if mu_hat>mu, otherwise 2(log L_hat-log L_mu).
One-sided asymptotic root q_mu=[Phi^-1(CL)]^2, CL>.5. Expected shape sensitivity
uses explicitly fractional background-only Asimov n=B. Sparse-count coverage is
unsupported/not validated; it is never presented as exact. Counting and shape
are different constructions, not interchangeable approximations.

A zero template produces infinity with explicit no_sensitivity status. Numeric
fits/roots are bracketed and checked. A nonlinear signal_counts(theta) uses an
explicit domain/grid and backend acceptance margin/metadata; all accepted
samples are retained, including disconnected regions. No arbitrary crossing is
reported as a unique upper limit. Between-point regions are unresolved.

The analytic NREFT convention and nuclear provenance documents are unchanged.
NREFT_CONVENTIONS.md SHA-256 remains
**a2c57658add45b273ca082d8b01fbab7dd739f8625d6f04a438d619d4cf385cd**.
No isotope masses/fractions, physical Ge responses, Wilson functions, nuclear
normalizations or physical detector/background values are added.

## Final validation and independent references

Runtime: Python **3.12.14**, NumPy **2.3.5**, SciPy **1.17.0**, pytest **9.1.1**.
Temporary pytest installation only; pyproject.toml is unchanged.

| Suite | Exact final result |
| --- | --- |
| Repository, `python -m pytest -q -W error` with temporary dependency path | **655 passed**, 0 failures, 0 skips, 0 warnings |
| Unchanged original pre-rate audit scripts | **154 passed**, 0 failures/skips/warnings |
| Unchanged rate-audit script, including original false-zero reproducer | **237 passed**, 0 failures/skips/warnings |
| Standalone synthetic example | successful; outputs counts, exact median CLs and approximate Asimov shape limits |
| Frozen document hash, dependency AST checks, git diff whitespace check | pass |

Final case accounting: 492 inherited +12 rate-product +62 detector +15
background +59 statistics +15 sensitivity = **655**. The three previously empty
analysis test files were filled; no inherited nonempty test/tolerance was altered.
External audit total is 391, separately 154+237.

New tests use elementary polynomial integrals, independent normal-CDF/PDF
references, independent finite Poisson sums/recurrences, and elementary likelihood
ratios. They do not obtain expected results from the production routine tested.
Tests cover identity/quenching Jacobian/count conservation, scalar/array/empty
energies, Lindhard zero and monotonicity, deterministic edges, narrow/broad/far-tail
Gaussians and partial-window losses, constant/varying truth and observed
acceptance, custom/matrix equivalence and ownership, binned exposure/edge rejection,
composite backgrounds, exact CLs and discrete medians, high/low backgrounds,
shape fits/ratios/one-sided thresholds, identical-shape splitting, high-count shape
improvement, explicit Asimov handling, impossible observations, invalid inputs,
root failures, warnings, immutable outputs, replaceable backends and nonlinear
scans. Import boundaries and all inherited missing nuclear-entry gates remain green.

The full pipeline oracle integrates a flat truth spectrum analytically against
a Gaussian via F(z)=z Phi(z)+phi(z), then applies explicit acceptance/exposure.
It derives B from flat rate×width×exposure, counting limits from independent
Poisson recurrences, and shape limits from elementary Asimov likelihood ratios.
A separate test adapts the existing synthetic NREFT fixture into a truth callable
and compares folded bins with independent SI constants and the closed polynomial
speed-PDF rate, without changing that rate kernel.

One new narrow-source test initially disagreed because the *fixture's* nominal
2*halfwidth differed from its actual rounded floating endpoints by about 8e-11
relative. Its density was corrected to 1/(upper-lower), explicitly normalizing
the supplied synthetic source. The analytic finite-width tolerance was not
relaxed. A narrow observed-window case is retained to verify the Gaussian truth
breakpoints prevent adaptive quadrature from missing its small support region.
No production clipping, warning suppression, xfail, skip, or tolerance weakening
was used to obtain the final results.

Reproduction commands:

```bash
# From the checkout:
PYTHONPATH=/tmp/wimp-legend-test-deps python -m pytest -q -W error
PYTHONPATH=src python examples/generic_analysis.py
# From the workspace (audit-repo symlink points to WIMP-LEGEND):
PYTHONPATH=/tmp/wimp-legend-test-deps:WIMP-LEGEND/src python -m pytest -q -W error audit-work/test_independent.py audit-work/test_falsifications.py audit-work/test_rate_audit.py
```

The supplied example uses synthetic r=2 on [0,5] keV_nr, q=.2, sigma=.1 keV,
truth/observed efficiencies .4/.6, flat b=3 and exposure 100 kg day. It returns
S approximately **(38.6291621332,95.7960516547,95.7963946196)** and
B **(60,120,180)**. Exact median counting mu_up approximately **0.1413434826**;
asymptotic Asimov shape mu_up approximately **0.1063046078**. These numbers
validate software assembly only; they are not physical sensitivity predictions.

## Files changed

- Rate closure: src/wimp_legend/rates/recoil.py,
  tests/test_rate_product.py, docs/RATE_AUDIT_002_RESOLUTION.md.
- Shared utility: src/wimp_legend/_analysis.py.
- Detector: all five existing detector/*.py files, tests/test_detector.py,
  docs/DETECTOR_RESPONSE.md.
- Background: backgrounds/__init__.py and model.py, tests/test_backgrounds.py,
  docs/BACKGROUNDS.md; radioactive.py/neutrinos.py unchanged.
- Statistics: statistics/__init__.py, counting.py, binned.py; new base.py;
  tests/test_statistics.py, docs/STATISTICS.md.
- Sensitivity: src/wimp_legend/sensitivity.py, new tests/test_sensitivity.py,
  docs/SENSITIVITY.md, examples/generic_analysis.py, this report.

All other tracked files are unchanged. No dependencies or physical datasets added.

## Remaining limits and independent-tester priorities

1. This is coder verification, not a new independent audit. Independently review
   the mathematical/statistical contracts, unit boundaries and numerical scope.
2. Exact counting assumes known background. Shape limits are asymptotic diagnostic
   results without sparse-count quantitative coverage. Collaboration likelihoods,
   nuisances and low-count toy calibration remain external/future work. Expected
   bands are omitted; expected shape is Asimov, counting is discrete median.
3. Adaptive integration needs caller-supplied narrow-feature/discontinuity points.
   Callable widths/custom models can introduce features the reference integration
   does not discover automatically. Validate convergence for intended analyses.
4. Matrices approximate sub-bin migration according to their declared averaging;
   no exact arbitrary-spectrum guarantee or implicit rebinning is made. Review
   efficiency ownership and truth coverage when supplying empirical models.
5. Nonlinear scans cover only explicit grid points. A unique upper limit or
   continuous parameter confidence set needs further assumptions/refinement;
   the default shape alternative fit is along each supplied point's shape.
6. Finite float support is not arbitrary precision. The final rate product fix
   cannot restore upstream erased values. Extreme quadrature/likelihood means,
   near-cancellation, subnormal precision and tiny-root tolerances remain bounded.
   Very large Poisson counts can be expensive or explicitly unsupported.
7. Pure callbacks/external objects may capture mutable caller state; immutable
   wrappers cannot certify their purity or freeze calibration versions. Record
   providers, configurations, source versions, exposures, ranges and numerical
   policies with analyses. Existing nuclear-source provenance remains external.
8. Physical Ge input, isotope/target counts, LEGEND calibration and backgrounds,
   and real-data validation all require separate sourced milestones. No normative
   nuclear convention or unresolved dataset normalization is resolved here.
