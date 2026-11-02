# Detector-folding / full-analysis speedup

## Delivery

Branch: `poster-speedup`, based on SI-speedup delivery
`0a4274686640a3e9e3d819091490a59f8916c9b8`.
Numerical implementation commit: `c2f1db68b4c9264a1e66cd9f690f243ccea736fb`;
tree: `8784199161ecce3a0591478e32d36684dca95403`.
The tree is identical to local verification checkpoint
`1c113414403bdd61944deab03f977e82dc1f5d9a`.
Subsequent report/evidence commits do not change numerical code.

## Implementation and unchanged science

* `folding.py`: optional bounded NumPy/SciPy Gaussian-bin evaluation, retaining
  signed CDF/survival differences, original batch-union pruning and original
  exact unpruned line probabilities. Arbitrary observed-efficiency callbacks
  keep the reference probability integral; custom response ports keep their
  own count implementation. No persistent dense response matrix or JIT.
* `models.py`: frozen `NoiseFanoSigma` replaces the existing noise-plus-Fano
  closure with the **same** formula and configured parameters. Only this known
  callback gets bulk evaluation; all arbitrary callbacks retain scalar semantics.
* `analysis.py`/`configuration.py`/master JSON: independent optional
  `numerical.folding_backend` (`reference` default, `vectorized` opt-in), selected
  for both signals and backgrounds. All quadrature and acceptance criteria are
  unchanged. Backend-specific fold keys prevent stale cross-backend cache reuse;
  truth-rate, exposure/style/on-grid-threshold reuse remains intact.
* Setup/main notebooks: opt-in example and display of resolved backend.
* Focused tests and `tools/benchmark_folding.py`: identical complete analysis
  inputs, independent cold caches/processes, warm reruns, saved numerical outputs,
  exact ROI checks, count/limit comparisons and per-process peak RSS.
* Prior `benchmark_si_fast.py`: migrate only the absent archived folding field
  to its reference default before its existing strict configuration equality gate.
  No physical parameter or assertion was removed from that gate.

There are **no changes under `src/wimp_legend/`**. Nuclear data/conventions,
SI/NREFT rates, quenching, resolution and efficiency formulas, background
physics/activities, detector/statistical kernels, ROI/CLs construction and
historical notebooks/data are unchanged. `poster`, `main` and `poster-workbench`
are not modified. No Numba, fastmath, extra workers or new dependencies.

## Profiling evidence

The complete two-mass O1 analysis profile (standard numerical settings, all
three scenarios/exposures and historical comparison) attributed 239.62 of
320.04 profiled seconds to the nested Gaussian-resolution probability routine;
270.98 seconds including its pruning wrapper. After the first vectorized
kernel, the same complete analysis took 76.96 profiled seconds and all saved
signals/backgrounds/24 limits were **exactly equal**. Scalar resolution sampling
still consumed 26.44 seconds; this motivated the final explicit array callable.
These are profiling measurements at the intermediate implementation
`30f8737d985c3746e4b462dc7df64bc9a8232ef6` (tree identical to local `8a4553b...`),
**not** the final unprofiled speedup claim.

A preliminary setup-overlapping profile is retained in evidence but excluded
from performance/physics acceptance; it loaded libraries while the environment
was being pinned. The repeated pinned baseline profile completed without stderr.

## Final unprofiled complete benchmark

| Complete 24-mass SI analysis | Reference | Vectorized |
|---|---:|---:|
| Cold wall time | 1533.126 s | 154.650 s |
| Warm wall time | 3.456 s | 3.434 s |
| Peak RSS | 214.72 MiB | 216.32 MiB |
| Warm cache misses | 0 | 0 |

Cold speedup: **9.9135x**. Warm speedup: 1.0064x (essentially unchanged). Both cold calls had 506 misses and two legitimate within-run hits; both warm calls had 508 hits.

All saved signals/backgrounds and **288 limits** compare exactly in this run: maximum relative differences are zero; ROI boundaries/statuses match and every original convergence check passes. No tolerance was relaxed.


Cold means a fresh subprocess with an empty numerical cache; OS page caches
are not flushed. Warm repeats exactly the same configuration in that process.
Both backends use `factorized_shm`, standard quadrature and convergence settings,
24 log-spaced masses 3–100 GeV, all five natural Ge isotopes, three scenarios,
all configured exposures and historical/revised comparisons. These are complete
SI analysis timings, not a raw-probability or raw-rate microbenchmark. Timed calls
exclude imports; subprocess time including imports is also recorded. Single
matched cold/warm samples on a shared 8-CPU-quota host; timings are not guarantees.
BLAS threads remain one. Peak RSS includes imports and both calls per backend.
There is no added JIT compilation overhead. The already-fast warm path may vary
in either direction because it reads caches rather than running CDF kernels.

## Validation

**37 focused tests passed in 176.58 s** on the final numerical tree.
The complete existing regression suite is **in progress**; no full-suite pass is claimed yet.
The earlier 35-test run is superseded/overlapping, not an additional test total.

No existing tests were modified, disabled or relaxed. New focused cases cover
all five isotopes at 3/10/100 GeV; folded counts and 90%-CLs limits; ER/NR continua,
lines, threshold bins and narrow/broader resolution; scalar/array/empty inputs;
constant and nonconstant efficiencies; arbitrary callback/custom-response/observed
fallbacks; zero-width rows; and ownership of response inputs/configuration.
Tolerances: probabilities exact in tested built-in cases; synthetic count
`rtol=2e-14`; physical count/limit `rtol=3e-8` (the earlier SI comparison tolerance).
The original coarse/fine convergence settings, significant-bin selection and
ROI-edge equality are retained.

## Limits and Tester action

See `FOLDING_SPEEDUP_TESTER_HANDOFF.md` and `docs/FOLDING_FAST.md` for commands
and independent checks. Independent acceptance remains pending. Warm caches
already avoid folding; remaining cold cost includes physical solar-source
integration, bin-crossing roots and output/cache I/O. No change to these physical
kernels or new interpolation of backgrounds was introduced to chase timings.
Existing provisional calibration/enrichment/background assumptions, source
closure deficits and physical nuclear-response validity uncertainties remain.
