# Background benchmark reports and tested source

The package implementation is the reviewed source from commit
`8c7261a213de7eae78bf711fdf29aa5e831fb3ff`, content tree
`a65512c12e8c4e48841dd5b193e42d394ebb2862`. This publication imports its
source bytes unchanged and retains the remote poster line-ending attributes.
GitHub publication through the connected app may use a different commit ID;
the IDs above identify the tested delivery, not a claim about publication
authorship or commit-object availability.

Start with [the latest independent report](CODE_AUDIT_008_BG_REMEDIATION.md),
then [the implementation response](CODE_AUDIT_007_IMPLEMENTATION_REPORT.md).
The earlier reports are preserved unchanged as historical evidence.

## Current acceptance

- Tester 008: PASS WITH CONDITIONS. The original Ar39/solar ES integration
  defects are closed; the tested public poster composite passes.
- Coder complete repository run: 889 passed, warnings treated as errors.
  The original log is included.
- The independently authored original suite: 183 passed unchanged, confirmed
  by Coder and Tester. The Coder re-run log is included.
- Tester repository run: 888 passed, one long-running check incomplete; this
  is distinct from the completed Coder run and is not added to its count.
- Tester new checks: 55 passed, 12 failed. Ten expose the inherited
  `limit-1` breakpoint workspace boundary; two expose writable private
  interpolation caches after deepcopy/pickle. These remain OPEN. No repair
  to these findings is included in this publication.

This is an illustrative benchmark, not a calibrated LEGEND prediction or
a quantitative neutrino-floor result. Physical/provenance limits and inherited
Ge limitations remain as described by the tester. Dense adaptive folds can be
slow; the plotting example uses the documented vectorized route.

## Work locally

```sh
python -m pip install -e '.[test,poster]'
python examples/background_benchmark.py --output bg-results
python -m pytest -q -W error
```

The full strict repository suite took about 28 minutes in the Coder environment.
See `docs/BG_BENCHMARK.md` for configuration, units and component replacement.
