# CODE_AUDIT_007 remediation — independent scientific acceptance

Review date: 2026-10-26. **Overall: PASS WITH CONDITIONS for the delivered numerical remediation and the tested background benchmark.** The four previously failing independent cases now pass unchanged. No production code or repository assertion was repaired or relaxed during this review. This acceptance does not establish a calibrated LEGEND sensitivity prediction.

| Area | Verdict | Finding |
| --- | --- | --- |
| Exact delivered source content | **PASS** | Both supplied patches independently reconstruct the declared tree. The claimed commit object itself was unavailable. |
| BG-AUDIT-001: Ar39 public counts | **PASS — CLOSED** | Full support/intensity and default Gaussian-window cases pass the original independent tests. |
| BG-AUDIT-002: solar ES public counts | **PASS — CLOSED** | Default Gaussian-window case passes against the unchanged independent raw-source reference. |
| Full public poster composite, including CEvNS | **PASS** | Actual default-policy API returns a finite, independently consistent count. |
| Shared integration change | **PASS WITH CONDITIONS** | Dense batching preserves features, warning behavior and global error accounting in tested cases. A pre-existing workspace-boundary exception remains. |
| Tabulated interpolation cache | **PASS WITH CONDITIONS** | Values and initial ownership are preserved; deepcopy/pickle lose the new private caches' read-only invariant. |
| Physical background model/provenance | **PASS WITH CONDITIONS, unchanged** | Prior source and physical-model limitations remain; no new physical input was supplied. |

## Exact source identity

Tested content tree: `a65512c12e8c4e48841dd5b193e42d394ebb2862`.

Coder's claimed commit: `8c7261a213de7eae78bf711fdf29aa5e831fb3ff`.

Incremental base: `5c62ee81923f75c4a753f8f6e0e8903eb183e7d1`, tree `48d3c6bd4fd55b558a6be0bd4d56f10637e9e794`. The incremental patch SHA256 is `eb2d0c27eb11f0b4b2d9045f6bff1aab5fdfe415b9f8a9bd76d421729d3748da`. Independently applying the cumulative patch (SHA256 `e45fcb0000479d9923f81288c8e9557f7173109a400dcbd3cb1c291c76b4a50f`) to `b2e6a423cb02acc64611d4e52dec047787345b85` produces the same tested tree.

The remote did not supply the claimed commit object, and direct SHA retrieval failed. Testing therefore used the patch-reconstructed index/work tree at the prior authenticated commit. `git write-tree` verifies the candidate content; `git diff --exit-code` verifies no unstaged deviation. This authenticates exact source bytes, not the unavailable commit's author, timestamp or parent metadata.

Only six files differ from the audited base: `_analysis.py`, `backgrounds/benchmark.py`, `backgrounds/beta.py`, two documentation files and the new `tests/test_bg_public_counts.py`. All 117 other tracked files, including the 22 existing test files, are byte-identical. Physical data tables, Ge responses, conventions, Wilson coefficients, rates, detector kernels, statistics and original notebook remain unchanged. Machine-readable identity and independent patch reconstruction are included in the evidence.

## Numerical and scientific comparisons

The source now supplies every Ar39 interpolation knot and the ES/CEvNS recoil endpoints of every fixed flux node. For dense breakpoint sets the helper integrates contiguous batches, retains each feature, allocates nonzero absolute tolerance by interval width, and sums results/errors with `math.fsum`. The combined error estimate is checked against the original global tolerance. This is a numerical dispatch correction; it does not renormalize a source or increase the requested tolerance.

The original independent suite remains byte-for-byte unchanged: test SHA256 `4bef22a188e4981c636b619195cc16df6d98082e0f781ff2ad4c9416dca543de`; reference-module SHA256 `79db64414cad6e2ce5cf8b1b6fde4294f24fb5bdfd084715522fd09dbbc6bbb7`. All **183 tests pass**, replacing 179 passes/four failures. Their references include an independent BetaShape run, direct piecewise-linear Ar39 integration, analytical Gaussian line probabilities, high-precision beta/decay equations, analytic ES integrals, independent solar kernels and absolute unit conversions. Production folding is not used to generate their expected results.

An additional actual public call uses the supplied poster configuration, Lindhard NR response, Gaussian/Fano resolution, efficiency 0.9 and the unchanged **default** policy (`epsrel=1e-8`, `epsabs=0`, `limit=300`). No plotting helper substitutes for this call:

| Quantity, 365.25 kg-day, 0.1–15 keV observed | Events |
| --- | ---: |
| `model.expected_counts([0.1,15],365.25)` | **83.08333162559286** |
| Prior independently integrated raw-source reference | **83.08333164312322** |
| Relative difference | **−2.10997 × 10⁻¹⁰** |

The public calculation takes about 137 seconds in this environment. Raw solar-spectrum integration and the production fixed-node representation are distinct approximations; agreement is not a claim of this many physically meaningful digits. The original independent ES reference is approximately 0.0384948281 events, versus 0.03849481047 for the fixed-node representation.

The regenerated example totals are 2918.49896017 events (natural Ge/AAr), 110.81592833 (enriched Ge/UAr), and 83.08333163 (poster). Poster Ar39, CEvNS and ES are respectively 0.01845753750, 0.41768132504 and 0.03849481047 events. The Coder's convergence tool completes; maximum component-integrated refinement change is about `1.04e-14`, and maximum reported significant-bin relative change is about `1.41e-11`. This is an internal quadrature check, separate from the independent raw-source comparison.

## Remaining numerical findings

**BG-RETEST-001 — inherited workspace-boundary limitation.** At exactly `limit−1` interior breakpoints the dispatcher still uses the old single `quad` call. Those prescribed intervals exhaust the workspace; QUADPACK can raise `IntegrationWarning` even for exactly integrable piecewise-linear data. Ten strict new cases fail at limits 2, 3, 7, 30 and 300, with both zero/nonzero absolute tolerance. All ten also fail unchanged against the audited baseline. This is an inherited explicit numerical exception, not a newly introduced wrong count or a reopening of the repaired dense Ar39/solar paths. The helper remains conditional on sufficient workspace for the requested feature set.

**BG-RETEST-002 — new minor private-cache invariant regression.** Fresh `TabulatedBetaShape` instances use read-only byte-backed arrays and resist write re-enabling. `copy.deepcopy` and pickle round trips instead produce writable owning arrays for `_energy_grid`/`_rate_grid`. The immutable public tuples and dataclass equality can then disagree with evaluation after direct private-cache mutation. With energies `(0,1,2)` and values `(1,2,0)`, changing the copied private midpoint to 200 changes evaluation at 1 from 0.8 to 80 while the raw tuple and equality with the original remain unchanged. Both strict round-trip tests fail. Copying alone preserves numerical values, and changing a private field is required to demonstrate divergent output; no ordinary poster count was found corrupted. This is nevertheless a real new failure of the advertised cache immutability contract. It remains open.

The new independent suite has **55 passes and 12 failures**: ten baseline workspace failures and the two new cache failures above. Passing cases cover dense-dispatch boundaries, exact trapezoid/antiderivative references on irregular grids, feature/support coverage, nonzero absolute-budget allocation, signed cancellation with a global acceptance/rejection check, warning propagation, nonfinite callbacks and scalar/vector interpolation equivalence. Expected integrals come from separate exact piecewise formulas and antiderivatives, not the production integrator. No failing assertion was removed.

## Regression evidence and acceptance limits

Repository regression status: **888 passed, one incomplete**. The unchanged suite with only `test_composite_counts_exposure_and_replacement` deselected completed with **888 passed, 1 deselected in 450.35 seconds** (`repository-other-final.log`). This includes all 14 new public-count regressions. A separate core run gives **830 passed** in 76.19 seconds; these overlap and must not be added to the 888. The full-suite attempt remained inside the older composite exposure/replacement check and was interrupted after more than 15 minutes without an assertion result (`repository-final.log`). It had completed 97 tests before that check. Consequently this review does **not** certify a completed 889/889 repository run or claim that the incomplete test passed or failed. The Coder's 889-pass claim is not independently confirmed here. Two earlier interrupted attempts also have incomplete logs; none establishes a numerical failure. Full regression completion remains an acceptance condition. The separate actual public poster-composite check and independent 183-case suite did complete successfully.

The unchanged pre-existing external background/SI/rate/falsification suites give **652 passes and one inherited failure** in 184 seconds. The failure is `test_dependency_boundaries`, which requires a single frozen-convention literal but finds the pre-existing duplicate in `wilson.py` and `germanium.py`. Neither file changed. The full original background independent suite gives **183 passes** in 376 seconds. The separate 67-case new suite and its ten-case baseline comparison are retained with their failures.

Environment: Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0, pytest 9.1.1, mpmath 1.4.1 and Matplotlib 3.10.8. Warnings are treated as errors. Timings are environment-specific.

The original physical/provenance conditions remain: only pp/B8 solar components; effective-density survival and atomic shell grouping; simplified target/transport and gamma/annihilation containment; assumed activities, 1% veto survival and sub-keV detector response; incomplete primary-PDF transcription verification, particularly weak Ga branches and Fe evaluation-version provenance. These inputs are unchanged, so this review reuses the original independent primary-source work rather than implying they were newly measured or fully verified. Prior Ge OBDME closure, diffraction-node roundoff and remote Gaussian-tail limitations are also unchanged. Ge-specific models were not introduced.

**Acceptance is scoped and conditional on completing the remaining repository check.** The original public Ar39 and ES defects are closed, and the tested poster public composite is numerically usable. Keep the workspace-boundary and serialization findings visible. Physical benchmark uncertainties and missing source/model input still prevent an unconditional quantitative LEGEND poster-sensitivity or neutrino-floor claim. The attached reproduction package contains strict tests, exact identities, comparisons and complete logs, including all retained failures.
