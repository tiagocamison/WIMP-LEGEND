# Folding speedup evidence

Numerical implementation: c2f1db68b4c9264a1e66cd9f690f243ccea736fb, tree 8784199161ecce3a0591478e32d36684dca95403.

`comparison.json` records the matched unprofiled 24-mass complete SI benchmark.
`reference/` and `vectorized/` each contain the saved run (manifest, spectra, limits, convergence) and timing metadata. Load either directory with `analyses.legend_poster.load_run(path)`; output checksums are preserved. Original absolute run paths in captured metadata are historical locations, not portable paths to execute.

`profiles/` contains full-analysis cumulative profiles of the intermediate implementation (local 8a4553b, tree identical to public 30f8737). Those timings include profiling overhead and are not the final speedup claim.

`execution/` retains bounded-run source/time/environment records and logs. Run 0001 overlapped dependency pinning and is excluded from performance/physics acceptance; run 0002 repeated the baseline in the pinned environment with empty stderr. Run 0003 is the initial overlapping 35-test pass; 0004 the complete intermediate profile; 0005 the final numerical-tree 37-test pass (local 1c11341, tree identical to public c2f1db6); 0006 the authoritative unprofiled complete benchmark on public c2f1db6.

Full regression is still running. No complete-suite pass is asserted in this checkpoint.

Reproduce with the commands in FOLDING_SPEEDUP_TESTER_HANDOFF.md. Sources, physical data, configurations and environment hashes are recorded in manifests; no cache files/pickle objects are needed for reading stored results.
