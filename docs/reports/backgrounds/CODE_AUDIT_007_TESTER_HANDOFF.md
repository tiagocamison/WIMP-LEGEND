# CODE_AUDIT_007 — independent re-review handoff

Candidate commit: `8c7261a213de7eae78bf711fdf29aa5e831fb3ff`

Candidate tree: `a65512c12e8c4e48841dd5b193e42d394ebb2862`

Use `identity.json` for the exact candidate commit/tree. The base is the audited
`5c62ee81923f75c4a753f8f6e0e8903eb183e7d1`.

## Restore

```sh
git clone --branch poster WIMP-LEGEND-poster.bundle WIMP-LEGEND
cd WIMP-LEGEND
python -m pip install -e '.[test,poster]'
python -m pytest -q -W error
python examples/background_benchmark.py --output bg-results
python tools/background_benchmark/check_convergence.py --output convergence.json
```

`CODE_AUDIT_007_FIX.diff` applies only to the audited commit. The cumulative
`BG_BENCHMARK_001_CORRECTED.diff` applies to the earlier
`b2e6a423cb02acc64611d4e52dec047787345b85` base. Apply only the appropriate
patch. Source tables retain their archival whitespace; use
`git apply --whitespace=nowarn PATCH`. Both exports must recreate the exact
candidate tree. A source ZIP and wheel are included for non-Git workflows.

## Re-run strict original failures

Use the poster configuration in `examples/bg_benchmark_config.json`, the supplied
Gaussian/Fano response and constant efficiency 0.9. Instantiate responses with
their **default** `IntegrationPolicy`: epsrel=1e-8, epsabs=0, limit=300.

1. Ar39: `component.expected_counts([-1,566],365.25)` with NoSmearing must
   reproduce mean source intensity times exposure/efficiency on full support.
2. Ar39: Gaussian observed `[0.1,15]` must agree with independent integration
   of the raw piecewise-linear table times Gaussian bin probabilities.
3. Solar ES: Gaussian observed `[0.1,15]` must return approximately
   0.03849481 events, consistent with the original independent raw-source
   reference 0.0384948281 within its documented quadrature representation error.
4. Full `model.expected_counts([0.1,15],365.25)` must return approximately
   83.083332 events. Re-test CEvNS separately: approximately 0.417681325 events.
   Its flux-node derivative kinks also required breakpoint metadata.

Use the original independent failing tests unchanged wherever available. This
handoff does not substitute for them or claim they have already been rerun.

## Original independent suite

The original evidence archive is included unchanged. Its historical README
refers to the old failing commit; set `BG_REPO` to the repaired checkout:

```sh
python -m pip install mpmath
unzip BG_BENCHMARK_001_INDEPENDENT_EVIDENCE.zip -d independent-evidence
BG_REPO="$PWD" PYTHONPATH="$PWD/src:$PWD/independent-evidence/audit7-work" \
  python -m pytest independent-evidence/audit7-work/test_background_benchmark_audit.py -q -W error
```

The archive can be unpacked in the checkout or adjust paths accordingly.
The unchanged suite now gives **183 passed**, replacing the original
179-pass/4-fail result. The Coder reran it; independent code re-review is still
requested. See `validation/independent-suite.log` and the original evidence
hashes in `identity.json`.

## Inspect numerical contracts

- Every Ar39 table knot and every ES/CEvNS flux-node recoil endpoint must be
  visible through `spectrum_breakpoints_keV`; shell/pp support markers remain.
- Large breakpoint sets must be sorted/deduplicated/filtered without dropping
  physical features, then partitioned into contiguous full-support batches.
- `limit` remains a per-call workspace bound. The same epsrel is supplied to
  each call; nonzero epsabs is allocated by width, not repeated per batch.
- The combined estimated absolute error must meet the original global target.
  Signed cancellation is tested; a combined overrun must raise, not return a
  count. Local IntegrationWarning and invalid callback values still raise.
- Frozen tabulated inputs and the new internal interpolation views must resist
  caller mutation and write re-enabling. Interpolation values are unchanged.
- No production source checksum, response domain, Wilson/NREFT convention,
  detector acceptance or nuclear normalization may change.

The complete repository suite gives **889 passed** (875 existing tests unchanged,
14 added). `tests/test_bg_public_counts.py` adds these checks without weakening old tests.
Analytic ES antiderivatives and beta table/CDF integration are independent of
the production count fold. Audit-number comparisons explicitly distinguish
source-representation accuracy from quadrature convergence.

## Scope and performance

Default scalar adaptive folding is deliberately conservative and can be slow
with thousands of truth knots and many observed bins. Test the public API
rather than using the faster plotting helper as a replacement. Acceptance
should be scoped to tested configurations, not arbitrary pathological detector
kernels or the inherited remote Gaussian-tail issue.

Retain the original audit's physical/provenance limitations and inherited Ge
limitations. In particular, these numerical corrections do not supply a complete
solar model, physical LAr transport, measured LEGEND activity/response, primary
PDF transcription completion or a quantitative sensitivity acceptance.
