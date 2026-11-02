# Fast SI validation evidence

`benchmark.json` records matched raw-rate timings, tracked allocation peaks,
35 pure-isotope folded comparisons and 84 SI-limit comparisons against the
unchanged generic standard snapshot. `fold_benchmark.json` records two cold
folds per algorithm with separate empty caches. Associated console logs are
included. The source and benchmark script hashes are recorded in the JSON.

The numerical benchmark source is `e76fd698a3f49d24679e4a88bc7662cb212494a9`.
The later iterator-ownership hardening at
`f93084236d7770a94e56eed911f47f04142d5e89` does not change numerical kernels or
the benchmarks' tuple inputs. The fold probe also records its uncommitted
validation-tool/test status at measurement time; executable physics hashes
identify the same numerical implementation. Do not conceal dirty status.

`full-suite.log` gives 1099 passes. It was collected before two extra regression
functions and the iterator ownership hardening. Post-hardening
`optimization-tests.log` gives 175 passes and `reference-tests.log` gives 219
passes. These overlap: **1101 distinct cases** are covered, not 1493, and not
a claim of one 1101-test invocation. `evidence.json` records this accounting,
source identities and evidence checksums. Earlier incomplete/failed development
attempts are not counted as passes.

Reproduce from the repository root in the documented editable environment:

```bash
python -m pytest -q -W error
python tools/benchmark_si_fast.py --output /tmp/wimp-si-benchmark
python tools/benchmark_si_fold.py --output /tmp/wimp-si-fold
```

The raw speedup excludes model construction and detector folding equally for
both algorithms. Cold-fold timings include them. Memory tracking is explicitly
scoped to Python/NumPy allocations; whole-process RSS is not a per-path native
memory comparison. Shared-host timings are not guarantees. The original
workbench physical assumptions and nuclear-data limitations remain unchanged.
