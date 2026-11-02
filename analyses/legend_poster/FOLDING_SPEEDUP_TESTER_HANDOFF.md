# Tester handoff — detector folding acceleration

## Candidate and authority

Task: accelerate the complete poster analysis on `poster-speedup`, retaining
accepted physics and the reference numerical backend. Independent acceptance
is pending; internal verification is not the Tester's acceptance decision.

Implementation commit: `12398ec1414fd76a40fd57b97f426eb84758ded2`.
Implementation tree: `126d4d04bf820669242c8f5de006332bc9dfa4a3`.
This tree exactly matches local checkpoint
`a8d90595fc4a38fc50e30f6e683f7b23aad00483`. Reports/evidence added subsequently
do not change the implementation. Base:
`0a4274686640a3e9e3d819091490a59f8916c9b8` (previous SI speedup delivery).

Only this branch is authorized for publication; no change to `poster`, `main`,
`poster-workbench`, physical datasets, frozen conventions or historical notebooks.

## Internal verification status

The complete 24-mass cold/warm benchmark passed, with exact saved-array/limit
agreement and a 9.9135x cold speedup. The 39 focused tests pass. The complete
regression suite passed **1140 tests in 2829.80 s**, warnings treated as errors.
The 39 focused cases overlap this total. Independent acceptance remains pending.

## Public workflow and interfaces

* Master JSON: `numerical.folding_backend` accepts `reference` (default; also
  absent-key fallback) and `vectorized`. It is independent of
  `signal.si_rate_engine`; enable `factorized_shm` as well for dense SI scans.
* `analyses.legend_poster.analysis.fold_component(component, bins, numerical,
  signal=False, refined=False)`: unchanged calling API. Backend selection is
  entirely numerical; events remain per kg-day at unit exposure.
* `analyses.legend_poster.folding.VectorizedGaussianBins(response, tail_sigma)`:
  same `bin_probabilities`/line interfaces as the reference wrapper.
  `prepare_folding(component, tail_sigma, backend='reference')` preserves routing.
* `run_analysis`, `load_run`, setup/main notebooks, cache and output formats are
  retained. The main notebook now displays the selected folding backend.
* `tools/benchmark_folding.py` runs separate fresh-process cold and warm complete
  analyses, independently checks outputs and writes timing/configuration/source
  records. A nonempty benchmark cache is an error, not a falsely labelled cold run.

## Mathematical checks the Tester should independently make

1. Gaussian bin integration: signed CDF differences for negative/crossing bins;
   survival differences for positive-tail bins. No window renormalization.
   NR mean uses the existing quenching, ER mean is the ER energy. Check the same
   observed-energy resolution evaluation point and shape/broadcast behavior.
2. Constant observed efficiency is applied exactly once. Truth efficiency remains
   outside the migration probability, multiplying the continuum/line source once.
   Nonconstant observed efficiency delegates the original probability integral.
   Only the known `NoiseFanoSigma` callback is called in bulk; arbitrary sigma
   callbacks still use core scalar sampling. Its formula and inputs are unchanged.
3. Every continuous spectrum uses exactly the same quadrature nodes, weights,
   interpolation/source breakpoints, bin crossings and pruning window. The window
   is the original batch-union window, not an additional per-event tail cut.
4. Discrete EC/other lines delegate the original exact, unpruned point-probability
   API. No continuum quadrature over delta functions. Observed spectra bypass
   both response and efficiency. Custom response ports keep their own count API.
5. Per-isotope counts, bulk composition weights, ROI edges and exact Poisson CLs
   must match. The SI rate selection gate, arbitrary signed proton/neutron
   coupling contraction and all nuclear source coefficients are unchanged.
6. Cache keys include selected backend and source/environment fingerprint. Changing
   only folding can reuse truth rates but not the other backend's folds. Changing
   exposure/style can reuse unit-exposure folds; new bin edges or detector noise
   must reintegrate. Corrupt cache arrays must still fail checksum checks.
7. Reproduce all five Ge isotopes at 3, 10 and 100 GeV and inspect thresholds/tails.
   Also inspect ER/NR continuum+line cases with noise 0.003/0.05 keV, nonconstant
   truth efficiency, scalar/multidimensional/empty inputs, zero-width rows and
   custom/nonconstant-efficiency fallback cases.

## Acceptance criteria and reproducible commands

Use Python 3.12 and `analyses/legend_poster/requirements.txt`; editable-install
this checkout. No Numba or new dependencies are required.

```bash
python -m pip install -r analyses/legend_poster/requirements.txt
python -m pip install -e .
python -m pytest -q -W error tests/test_folding_fast.py
python -m pytest -q -W error
python tools/benchmark_folding.py --output /tmp/fold-profile --masses 2 --profile-analysis
python tools/benchmark_folding.py --output /tmp/fold-dense --masses 24
```

Use fresh output directories; existing caches are not silently emptied.
The profile run is for cost attribution, not the unprofiled speedup claim.
Dense scan uses standard quadrature settings, 24 log-spaced masses from 3 to
100 GeV, all three configured scenarios, all five Ge isotopes, all exposure
settings and the historical comparison. Both backends use factorized SI rates.

Required comparisons:

* Gaussian array probabilities: exact array equality against the reference in
  the focused built-in response cases; weighted-efficiency fallback equality.
* Synthetic continuum+line folds: `rtol=2e-14`, `atol=peak_bin*1e-15`.
* Physical isotope folds and benchmark outputs: `rtol=3e-8`,
  `atol=peak_reference_bin*1e-13`. Limits and mu: `rtol=3e-8`, `atol=0`.
* Exact ROI-edge equality, same result status and all original coarse/fine
  convergence gates. Existing tolerances and tests must not be weakened.
* Review cold/warm wall times separately. Per-process peak RSS includes imports
  and both calls; compilation overhead is zero. No oversubscribing workers,
  fastmath or persistent response matrices are introduced.

## Evidence and remaining assumptions

See `FOLDING_SPEEDUP_REPORT.md` and `validation/folding_speedup/` for results,
logs, recorded configurations, source hashes and timing evidence. Reports must
state actual execution counts; overlapping targeted/full tests are not additive.
The existing benchmark assumptions remain: provisional enrichment/calibration
and backgrounds, AFH closure deficits, unvalidated physical high-q range and
known-source statistical assumptions. This optimization does not resolve them.

For a separate independent pass, concentrate on fallback selection, tail
cancellation, deterministic bin boundaries, cache invalidation and complete
limits rather than accepting only a fast raw probability microbenchmark.

## Zero-width hardening / validation in progress

The corrected implementation delegates any zero-width Gaussian batch to the
reference wrapper, preserving its local final-bin endpoint semantics. Two new
ER/NR edge regressions pass. This does not affect the recorded poster benchmark,
whose resolution is strictly positive; its measured source was C2 (`c2f1db6`).
The full suite was deliberately stopped on that earlier revision after this
review finding. The cancelled run (exit 130) is preserved under
`validation/folding_speedup/execution/run-0007/`; it is not counted as a pass.
The full suite passed on `12398ec`: **1140 passed in 2829.80 s**.
The matched complete benchmark finished; see the final delivery update below.
Independent acceptance is pending.

The completed full-suite record is preserved under
`validation/folding_speedup/execution/run-0010/`. Run 0009 lost its executor;
its partial log is preserved separately and is not counted as a pass.

## Final delivery update

The corrected source above passed **1140 tests in 2829.80 s**, including all 39
focused tests. The final 24-mass repeat measured **9.8954x cold speedup**, with
exact agreement of saved arrays and all 288 limits. Warm runs were 3.271/3.393 s;
peak RSS was 208.59/208.27 MiB; added JIT compilation was zero.

The final raw benchmark/profile records did not survive the interrupted upload's
workspace rollback. `validation/folding_speedup/FINAL_MEASUREMENTS.json` is an
explicitly recovered tool-completion summary, not a fabricated raw manifest.
The final full-suite log and earlier complete benchmark raw data are durable.
Reproduce the final-source paired benchmark with the commands above if raw
final-source timing evidence is required for independent acceptance. No additional
benchmark runs are required to pull/use the fully updated implementation.
