# Detector folding speedup — final report

Numerical implementation: `12398ec1414fd76a40fd57b97f426eb84758ded2`.  
Numerical tree: `126d4d04bf820669242c8f5de006332bc9dfa4a3`.
Branch: `poster-speedup`. All implementation code is pushed.

## Result

**1140 tests passed in 2829.80 seconds**, with warnings treated as errors.
This single complete invocation includes all **39 focused folding tests**;
overlapping test counts are not added. The full log and immutable source record
are published in `validation/folding_speedup/execution/run-0010/`.

The final 24-mass complete SI analysis measured:

| Measurement | Reference | Vectorized |
|---|---:|---:|
| Cold wall time | 1568.683 s | 158.526 s |
| Warm wall time | 3.271 s | 3.393 s |
| Peak process RSS | 208.59 MiB | 208.27 MiB |
| Warm cache misses | 0 | 0 |
| Added JIT compilation | none | none |

**9.8954x faster cold.** All saved signal/background arrays and **288 limits**
matched exactly; maximum relative differences were zero. ROI boundaries,
statuses and existing convergence criteria matched. Warm performance is
essentially unchanged because those runs read numerical caches.

Both modes used factorized SI, standard numerical settings, 24 logarithmically
spaced masses from 3 to 100 GeV, all five Ge isotopes, three scenarios, configured
exposures and historical comparisons. Separate processes and empty numerical
caches defined cold runs; OS caches were not flushed. BLAS threads were one.
RSS includes imports and both calls. These are shared-host measurements, not
runtime guarantees.

## Implementation

* Optional NumPy/SciPy Gaussian-bin batches replace nested Python loops, retaining
  the reference CDF/survival subtraction and existing integration nodes/pruning.
* Immutable `NoiseFanoSigma` evaluates the existing noise-plus-Fano formula in
  bulk. Arbitrary sigma callbacks retain scalar sampling; nonconstant observed
  efficiencies and zero-width Gaussian batches delegate reference probabilities.
* Discrete lines retain exact unpruned probabilities; observed templates bypass
  response/efficiency. Custom response ports retain their own count interface.
* Backend-specific fold keys preserve cache integrity; existing isotope and
  unit-exposure reuse remains. No Numba, fastmath or worker pool was added.
* Master configuration, notebook examples, regression tests and reproducible
  `tools/benchmark_folding.py` are included. Package physics modules, datasets,
  conventions, statistical methods and historical notebooks are unchanged.

Profiling identified nested Gaussian probabilities as the dominant removable
cost. The final optimized two-mass complete profile took 50.205 profiled seconds;
profiling overhead is separate from the dense timings above. Existing tests and
acceptance tolerances were not relaxed. Two ER/NR local-edge regressions cover
the zero-width fallback; interrupted attempts are retained separately from passes.

## Enable and reproduce

Set in the saved master configuration:

```python
config['signal']['si_rate_engine'] = 'factorized_shm'
config['numerical']['folding_backend'] = 'vectorized'
```

Use the setup/main notebooks normally; restart kernels after pulling code.
The reference implementations remain defaults. Reproduce the full benchmark:

```bash
python tools/benchmark_folding.py --output /tmp/folding-benchmark --masses 24
```

Use a fresh output directory. See `docs/FOLDING_FAST.md` and
`FOLDING_SPEEDUP_TESTER_HANDOFF.md` for APIs, equations and independent checks.

## Evidence availability

The complete test log and earlier matched benchmark raw numerical data are
published. The final benchmark/profile finished successfully, but their raw
local records were lost after the interrupted upload reverted the workspace
snapshot. `validation/folding_speedup/FINAL_MEASUREMENTS.json` faithfully records
the observed completion measurements and this limitation; it is not a substitute
raw run manifest. No additional calculations were performed for this resupply.
Independent Tester acceptance remains pending. Existing provisional physical
background/calibration/enrichment and AFH validity assumptions are retained.
