# CODE_AUDIT_007 — implementation response

Candidate commit: `8c7261a213de7eae78bf711fdf29aa5e831fb3ff`

Candidate tree: `a65512c12e8c4e48841dd5b193e42d394ebb2862`

The independent audit is accepted: the original poster plotting calculation
passed conditionally, but the normal default public count API failed. The
original 875 passing tests did not establish correctness of that path, and the
statement that only physical uncertainties remained was too strong.

## Corrections

- Ar39 now supplies all 5,651 interpolation knots to the existing truth-response
  fold. It remains the same normalized pinned BetaShape table.
- Solar ES now supplies every fixed flux-node recoil endpoint, plus the existing
  electron-shell edges. A new default-composite regression exposed the same
  convergence issue in CEvNS; its node endpoints are now supplied too.
- The shared `IntegrationPolicy.integrate` handles more breakpoints than a single
  QUADPACK workspace by integrating contiguous bounded batches. Relative
  tolerance, full support, warning-as-error behavior and per-call workspace
  limit remain unchanged. Absolute error allowance is shared by interval width.
  The sum of estimated errors must satisfy the original global tolerance;
  otherwise an `IntegrationWarning` still raises. There is no fallback that
  accepts a failed integral or discards a source feature.
- Immutable-byte-backed interpolation views avoid repeatedly converting the
  entire Ar39 table for each scalar quadrature call. They cannot be made
  writable and do not alias caller input. Public raw tuples and equality remain
  unchanged.

No public API signature changed. Detector, rates, statistics, interactions,
Ge nuclear data, background input snapshots, frozen NREFT document and the
historical notebook remain unchanged. This is a numerical integration repair,
not a physical-model update or normalization change.

## Verification

Before correction, the three audit cases reproduced: Ar39 full support with
NoSmearing, Ar39 Gaussian/Fano poster window, and solar ES Gaussian/Fano poster
window all raised IntegrationWarning. The new tests use default policies;
previous tests and their tolerances are untouched.

Fourteen regressions cover dense-breakpoint integrals, workspace/absolute-budget
allocation, combined-error rejection including signed cancellation, warning
propagation, buffer immutability, full-support Ar intensity, Gaussian Ar bins
against direct table/CDF quadrature, unbroadened ES against analytic polynomial
integrals, and Gaussian ES/CEvNS plus the entire composite against the supplied
independent audit. They do not use the plotting helper to produce expectations.
The audit comparisons allow the documented raw-flux versus fixed-quadrature
representation difference; the integrator's tolerance is never relaxed.

The auditor's original `BG_BENCHMARK_001_INDEPENDENT_EVIDENCE.zip` was also
recovered and hash-verified. Its scientific test file was rerun unchanged:
**183 passed in 327.27 s**, including all four original failing cases. This is
a Coder-run reproduction of independently authored tests, not a new independent
review verdict. Original source/test evidence is included without modification.

The complete repository suite passed: **889 passed in 1658.22 s (27m38s)**
with warnings treated as errors: 875 existing tests unchanged, 14 added.
Exact commit/tree identity and the test result are in `identity.json` and
`validation/final-full-suite.log`. The wheel passes loading/dense-integration
smoke checks with both background and Ge source snapshots packaged.

Regenerated poster figures use the corrected source breakpoint metadata. For
365.25 kg-day in 0.1–15 keV observed, the total remains approximately 83.083332
events. Agreement is numerical agreement of this benchmark, not physical
predictive precision. Updated convergence results and component counts are
included. Immutable source files and the background manifest match the audited
commit byte for byte; no source regeneration was needed for this repair.

Dense adaptive public folding is slower than the vectorized plotting path,
especially for many observed bins. That remains a practical performance
limitation, rather than a reason to loosen accuracy or bypass the count API.
Arbitrary user kernels/responses can still require their own breakpoints; the
existing detector's remote Gaussian-tail limitation was not repaired here.

## Remaining acceptance conditions

Independent re-review is required to close the audit's public-count FAIL; this
Coder report is not an independent acceptance. Physical conditions are unchanged:
pp/B8-only solar coverage, approximate target/atomic response, activity history,
source-shaped LAr transfer, assumed veto/containment, and uncalibrated low-energy
response. The incomplete primary-PDF branch transcription review and evaluation
version gaps remain explicit. No quantitative LEGEND prediction or neutrino-floor
claim is granted. Inherited convention-literal duplication, OBDME closure and Ge
diffraction-node negativity are outside this bounded fix.
