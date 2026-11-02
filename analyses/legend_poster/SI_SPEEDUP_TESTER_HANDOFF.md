# Tester handoff — optional fast SI rate

## Revision and current status

Branch: **`poster-speedup`**, based on completed workbench delivery
`f5863ad280e0fc2250b6926d1659920e4f1caf5a`.

Numerical implementation and timing-source commit:
**`e76fd698a3f49d24679e4a88bc7662cb212494a9`**.
Numerical implementation tree: **`b8315e0a08bd388c19e3e9190ce543b201d79a3b`**.
This matches locally tested implementation `890eeffeb12a9a6de0d7e06b44272295906fa268`.
Subsequent validation commits add this document, the cold-fold driver and
regressions. Input-ownership hardening is published at
`f93084236d7770a94e56eed911f47f04142d5e89`, tree
`d499813d4a577651dc4d390acdbb00fc66551426`: breakpoint iterators are copied before
matching the explicit SHM. The numerical kernels and benchmark tuple inputs are
unchanged.

Status on 2026-10-26:

* **175 focused optimization tests passed in 154.44 s**. This includes the
  additional folded-count/90%-limit test for four coupling cases and iterator ownership.
* **219 existing SI/rate/halo tests passed in 38.45 s**; these overlap the full
  suite and must not be added to its eventual total.
* **35 pure-isotope folded spectra and 84 SI limits** passed comparison with
  the immutable generic standard snapshot. Maximum tested-bin relative
  difference in the first benchmark was `1.9306778398231472e-13`; maximum
  relative limit difference was `4.440892098500626e-16`.
* The matched cold detector-fold probe completed: median speedup **3.2586x**
  for the configured middle mass (10 GeV), Ge-73, standard numerical settings.
  The final raw-rate scan measured **2116.45x**, not an end-to-end speedup.
  Both samples, allocation measurements and source hashes are included in
  `validation/si_speedup/`.
* The **complete regression run passed: 1099 tests in 4537.31 s**, with warnings
  treated as errors. It was collected before two extra regression functions and
  the iterator input hardening. Post-hardening 175 optimization and 219 existing
  affected-interface tests passed in separate runs. They overlap: **1101 distinct
  cases** are covered, not a single 1101-test invocation. No existing test or
  assertion was edited or disabled. The final checkout should collect 1101.


This is an implementation-author handoff, **not independent Tester approval**.
The final implementation report and checksummed evidence record these completed runs.

## Checkout and reproduce

```bash
git fetch origin
git switch --track origin/poster-speedup   # omit --track if it already exists
python -m venv .venv
.venv/bin/python -m pip install -r analyses/legend_poster/requirements.txt
.venv/bin/python -m pip install -e '.[test,workbench]'
.venv/bin/python -m pytest tests/test_si_fast.py -q -W error
.venv/bin/python -m pytest -q -W error
.venv/bin/python tools/benchmark_si_fast.py --output /tmp/wimp-si-benchmark
.venv/bin/python tools/benchmark_si_fold.py --output /tmp/wimp-si-fold
```

Run from the repository root. The complete suite contains a slow nested
background-count test; retain its final result rather than terminating it and
claiming the partial log passed. Performance measurements take minutes and
depend on host contention. The benchmark requires the default physical master
configuration: it fails rather than comparing a changed physical configuration
to a frozen reference as though they were identical.

Enabling the optimization in `analyses/legend_poster/config.json`:

```python
working['signal']['si_rate_engine'] = 'factorized_shm'
```

Validate/save with `poster_setup.ipynb`, then rerun `poster_workbench.ipynb` or
`python -m analyses.legend_poster --profile standard`. `"generic"` is the
default; an old config without the key also stays generic. Both active notebooks
expose the setting. Keep the original standard snapshot unchanged.

## Interfaces to audit

* `wimp_legend.halo.shm_eta.mean_inverse_speed_shm(halo, v_min_c)`.
* `wimp_legend.si_fast.is_o1_only(coefficients)`.
* `wimp_legend.si_fast.differential_rate_shm_per_kg_day_keV(E_nr_keV,
  halo=..., use_fast_si=False, **generic_arguments)`.
* Existing `SIO1Signal`, extended only with optional `use_fast_si=False` and
  `halo=None`. A fast request must supply the matching SHM PDF, support and
  breakpoints; contradictory declarations raise.
* Workbench `signal_model`, `pure_signal_counts`, and the run-local SI limit
  adapter reuse in `analyses/legend_poster/analysis.py`.
* `tools/benchmark_si_fast.py`, `tools/benchmark_si_fold.py`, and
  `tests/test_si_fast.py`.

The generic rate engine, particle responses, frozen convention document,
accepted AFH datasets, target handling and detector/statistical kernels are
unchanged. No Helm response or new physical nuclear/background data is present.
Detailed derivation/API notes are in [`docs/SI_FAST.md`](../../docs/SI_FAST.md).

## Mathematical oracle and units

Starting from the existing generic assembly, independently derive

\[
C_M(q^2)=\sum_{\tau,\tau'=0,1}c_1^\tau c_1^{\tau'}W_M^{\tau\tau'}(q^2),
\qquad
\frac{dR}{dE_{\rm nr}}=
\frac{N_T\rho_\chi}{m_\chi}\,c\,(86400)\,
\mathrm{GEV\_MINUS2\_TO\_CM2}\,(10^{-6})
\frac{2m_T}{2J+1}C_M(q^2)\eta(v_{\min}).
\]

Use the frozen half-isospin convention: `c0=(cp+cn)/2`, `c1=(cp-cn)/2`.
Masses are GeV, Wilson coefficients GeV^-2, q^2 GeV^2, recoil input keV_nr,
speeds beta=v/c and eta an integral against d(beta). Output is
events/(kg day keV_nr). The four ordered M entries must all contribute, with
their signs intact. Do not transpose, symmetrize, take absolute values, or
repair source particle-number closure deficits.

For the halo, compare to direct integration of the accepted boosted/truncated
speed PDF. Freese/Lisanti/Savage arXiv:1209.3339v3 Appendix B Eq. (B4), hard
truncation parameter beta=0, is the analytical reference. That paper's truncation
parameter is distinct from our speed variable. Narrow exponential-difference
integrals use 16-node Gaussian evaluation to avoid subtraction cancellation;
there is no unboosted approximation or hidden cutoff.

For direct physical SI normalization, independently check
`sigma_p = GEV_MINUS2_TO_CM2 * mu_chi_p**2/pi * cp**2` and
`sigma_p_90 = mu_90 * sigma_p_ref`. The supplied numerical proton/nucleon
reference masses and all conventions remain explicit. Pure-neutron models use
the coefficient dispatcher; they cannot generally be encoded by a finite cn/cp
ratio and zero reference proton cross section.

## Required independent checks

1. **Selection gate:** with opt-in enabled, test every non-O1 operator and both
   isospin indices, including the smallest nonzero float. Verify the generic
   function is called with the entire coefficient object intact. A nonzero
   coefficient must not be ignored because its square underflows, the target
   has zero spin, j_chi=0, the recoil is inaccessible, or a response happens to
   vanish. With opt-in disabled, even pure O1 must call generic. All-zero
   coefficients qualify for the fast path and give zero rates.
2. **Raw physics:** all stable Ge isotopes 70/72/73/74/76, at least masses
   3/10/100 GeV and E=0, interior, near-support and outside-support recoils.
   Test cp=cn, unequal cp/cn, opposite signs, global sign reversal, proton-only,
   neutron-only, both zero and ratios close to destructive cancellation.
   Use the *existing generic integrator* as oracle, not a copy of the fast code.
3. **Ordered terms:** include a clearly synthetic fixture with different 01
   and 10 entries to detect accidental doubling/symmetrization. Retain signed
   contractions and check the nuclear spin denominator and all unit factors.
   Independently verify the free-proton/coherent normalization if possible.
4. **Eta stability:** normal, zero/tiny boost, boost at/above escape speed,
   small escape sphere, negative/infinite/NaN thresholds, scalar/array/empty
   shapes, monotonicity and exact zero at/above support. In the last tiny tail,
   adaptive double-precision relative-only quadrature is not a reliable oracle:
   use high-precision angular/PDF integration. A rounded E_max can map one ulp
   below vmax; require the same accessibility predicate as generic, not an
   artificial energy clamp. Check that no source warning is suppressed.
5. **Detector/statistics:** reproduce the 35 isotope folds (seven configured
   masses times five isotopes) and all 84 saved O1 limits, including natural/
   provisional-enriched targets, exposures and historical/revised backgrounds.
   Repeat folded comparisons with isospin violation, pure neutron/proton and
   zero couplings, using synthetic background counts to isolate the rate change.
   Check ROI edges, retained fractions, event counts and the exact 90% Poisson-CLs
   construction. Confirm neither efficiency nor detector response is repeated.
6. **Caches:** exposure, background activities/settings and a threshold moved
   to an existing bin edge must reuse pure/folded signals without invoking the
   truth-rate callback. Changing NR noise requires new folds but can reuse
   truth rates. An off-grid threshold changes actual bin edges and therefore
   legitimately requires new integrals. Switching the SI option must not
   invalidate non-SI results. Composition changes must use the existing correct
   per-kg isotope weights. Source/config/data changes must not reuse stale data.
7. **Validation gates:** incompatible nuclear convention/units, missing required
   entries, invalid energies/spins/masses and a contradictory SHM declaration
   must fail. Inspect active M domain/validated-range enforcement. Generic still
   evaluates all its channels; factorization evaluates the relevant M entries.
   Do not silently broaden any source-supported range.
8. **Performance:** benchmark identical vectors, model inputs and precision.
   Separate raw rate timing from complete cold detector folds and cache hits.
   The raw benchmark is 12 masses x 5 isotopes x 24 recoil energies, twice,
   with model construction excluded equally from both paths. Cold-fold timing
   includes construction, truth-grid refinement and detector integration with
   separate empty caches. Record both samples, not only the most favorable one.

## Numerical acceptance budgets

These are implementation-comparison budgets, not physical uncertainty estimates:

| Quantity | Relative tolerance | Absolute tolerance |
| --- | --- | --- |
| Eta versus adaptive PDF moment, away from an unresolvable terminal interval | 2e-9 | 2e-14 times peak reference eta |
| High-precision terminal-tail oracle at fractional distances 1e-4 and 1e-7 | 1e-7 | 1e-25 |
| Raw AFH rates | 3e-9 | 2e-13 times peak absolute reference rate in the test vector |
| Saved-standard folded spectra | 3e-8 | 1e-13 times peak reference bin count |
| Extra signed/neutron-only folded cases | 3e-9 | 2e-13 times peak reference bin count |
| Physical SI limits against saved standard | 3e-8 | 0 |
| Extra coupled-case mu_90 | 3e-9 | 0 |

The terminal-tail budget accounts for float64 support endpoints; its absolute
contribution is negligible, but must remain finite, nonnegative and supported.
Keep the standard workbench's separate interpolation/quadrature convergence
checks too. Do not relax its existing tests to make this optimization pass.

Memory evidence uses `tracemalloc` per 24-point call (tracked Python/NumPy
allocations, excluding untracked native memory). Whole-process peak RSS includes
both algorithms and imports and is not a per-algorithm native-RAM comparison.
Do not describe either metric as total memory saved without further measurement.

## Remaining assumptions and evidence delivery

All previous target/calibration/activity/halo/nuclear-range uncertainties remain.
The fast path is for the explicit SHM and O1 only; arbitrary PDFs and other
operators retain generic integration. Detector integration is still a substantial
cold-scan cost. This change does not validate LEGEND calibration, nuclear inputs
or earlier independently documented core audit limitations.

Final benchmark JSON and regression logs are in `validation/si_speedup/`; the
implementation report is `SI_SPEEDUP_REPORT.md`. Check their source commits,
actual file/data hashes, reference snapshot hash and resolved configuration.
The evidence records collection timing and overlapping test counts explicitly.
Record independent failures and tolerances rather than replacing the supplied
assertions with copies of the implementation.
