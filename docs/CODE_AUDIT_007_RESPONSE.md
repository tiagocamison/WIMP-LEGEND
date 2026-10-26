# CODE_AUDIT_007 — background public-count correction

Review input: `CODE_AUDIT_007_BG_BENCHMARK.md` (2026-10-25), reviewing
commit `5c62ee81923f75c4a753f8f6e0e8903eb183e7d1`, tree
`48d3c6bd4fd55b558a6be0bd4d56f10637e9e794`.

## Findings and bounded fixes

**BG-AUDIT-001 accepted.** Before correction, Ar39 public counts with the default
`IntegrationPolicy` raise `IntegrationWarning` both over full support with
`NoSmearing` and for the poster `[0.1,15]` observed-keV Gaussian/Fano window.
The spectrum is piecewise linear, but its 5,651 source knots were not supplied
to the response. `build_benchmark` now passes those exact knots through the
existing `TruthBackground.spectrum_breakpoints_keV` interface. It does not
smooth, coarsen or renormalize the source table.

**BG-AUDIT-002 accepted.** The fixed solar flux quadrature creates a distinct
ES recoil endpoint for each node, in addition to atomic-shell jumps. Passing
only each source's highest endpoint left those kinks hidden. The factory now
passes every ES node's exact finite-mass recoil endpoint, together with the
existing shell and narrow-support breakpoints. The new default composite test
also reproduced the same kind of convergence failure for CEvNS derivative
kinks. All CEvNS quadrature-node endpoints are now supplied as well. The
scattering formula, flux weights, probabilities and atomic-shell approximation
are unchanged.

Both sets exceed one QUADPACK workspace. `IntegrationPolicy.integrate` now
partitions such a set into contiguous batches, each using at most half the
workspace for prescribed intervals and reserving the rest for adaptive
refinement. It retains `epsrel`, `limit`, warning-as-error behavior and the
full integration domain. A nonzero absolute error budget is divided by batch
width; it is not repeated in full for every batch. Values are combined with
`math.fsum`, and their combined estimated absolute error is checked against
the original global tolerance; an
excess still raises `IntegrationWarning`, including for signed cancellation.
Smaller breakpoint sets continue through the original single-call path. No exception is caught and
accepted as a numerical result; no tolerance is relaxed.

This shared numerical helper required a small implementation change to support
the source-provided dense breakpoint sets. Public signatures, detector code,
rate/statistical engines, nuclear-response data, conventions and the historical
notebook are unchanged. No new response or plotting-only count API is introduced.

The tabulated beta implementation also retains two immutable-byte-backed NumPy
views of its owned input tuples, avoiding full-table reconversion at each
scalar quadrature call. Interpolation and source values are identical. Views
cannot have write access re-enabled, and do not alias caller-owned arrays.
There is no change to production source files or their hashes.

## Regression gates

`tests/test_bg_public_counts.py` adds default-policy gates without modifying
previous tests or their explicitly loose policies:

- Dense piecewise-linear integration at workspace limits 1, 2, 7 and 300,
  including unordered, repeated and out-of-range breakpoint handling, plus
  rejection of a combined error estimate exceeding the global target.
- Instrumented quadrature verifies unchanged relative tolerance/workspace,
  contiguous domain coverage and a shared absolute error budget.
- Integration warnings still raise; immutable interpolation buffers remain
  protected and the tabulated object's equality semantics are preserved.
- Full-support Ar39 with no smearing reproduces independently computed
  run-averaged source intensity times exposure and efficiency.
- Gaussian Ar39 bin counts agree with direct Gauss-Legendre integration of
  each raw linear table interval against independently evaluated Gaussian CDFs.
- Unbroadened ES agrees with analytic integration of each node's polynomial
  differential cross section above each electron binding threshold.
- Gaussian ES, CEvNS and the complete public `BenchmarkModel.expected_counts` call
  agree with the independent audit's poster-window results. Those comparisons
  allow the audit's documented raw-flux versus fixed-quadrature representation
  difference; they are not a relaxation of the production integrator tolerance.

No regression reference calls the plotting fold. The exact run results and
source/commit identities are provided in the follow-up delivery report.

## Acceptance limits retained

The earlier statement that only physical uncertainties remained was too strong:
the original 875 passing tests did not establish a working default public fold.
The public-count regressions now test the failing path explicitly. These fixes
require independent re-review before claiming the audit's public-count FAIL
has been independently closed.

The original audit's physical conditions remain: pp/B8-only solar coverage,
benchmark activities/history, source-shaped LAr-to-Ge transfer, assumed veto,
all-or-nothing containment, grouped L1/M1 energies, effective target surrogates
and uncalibrated low-energy response. Full primary-PDF transcription/provenance
review is incomplete, especially minor Ga68 branches and Fe55 evaluation dating.
No coefficient or physical input was added or altered to address these gaps.

Inherited convention-literal duplication, Ge OBDME closure and Ge70/72
float64 diffraction-node negativity are outside this correction. No NREFT
normalization, data rescaling or nuclear-response repair is attempted. This
remains an illustrative benchmark, not a validated LEGEND prediction or a
quantitative neutrino-floor calculation.

Dense adaptive public folds are a correctness/reference path and can be slow
for many observed bins; this correction does not promise interactive plotting
performance. The existing vectorized plotting example remains available.
Acceptance is scoped to the tested responses/windows: arbitrary user-supplied
kernels or pathological efficiency/resolution functions can still require their
own breakpoints or produce explicit integration errors. The detector's inherited
remote Gaussian-tail limitation is not re-audited or repaired here.
