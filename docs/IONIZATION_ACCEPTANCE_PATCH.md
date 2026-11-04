# Ionization acceptance patch application

Branch: `feature/pre-smearing-acceptance`.
Numerical code commit: `f82f84c4b4d3662cd1536143e1a7a1aa4cac86f6`.
Tested tree: `f380defcbcb09edb158b893505250f44a89fb1af` (identical to local test checkpoint `5e744a2`).

The supplied patch had stale context: an extra existing Gaussian-breakpoint
comment and terminal blank-line differences prevented direct `git apply`.
Its intended edits were applied while retaining the existing comment and plots.

`detector.presets.<preset>.ionization_acceptance` supports `none`, `hard`, and
`erf`. Acceptance is evaluated before smearing, on ER truth energy directly or
on the existing quenched mean NR visible energy. Discontinuities/narrow turn-ons
are supplied to both quadrature paths in truth-energy units. Rejected continua
and lines never invoke the smearing kernel. Already-observed templates retain
their direct route. Default `none` preserves the previous selection behavior.
The supplied 0.020/0.005-keVee settings are disabled placeholders, not measured
LEGEND calibration data. Choose physical parameters explicitly before enabling.

Added direct-construction validation and regressions for narrow turn-ons,
pre-smearing rejection and input validation. The branch's missing explicit
`numerical.folding_backend` was restored to `reference`; this matches its prior
implicit default and fixes an existing test without changing runtime selection.

Validation: **121 passed, 15 deselected in 63.01 seconds**, warnings as errors.
All 19 ionization tests passed. Command:

```bash
python -m pytest -q -W error tests/test_ionization_acceptance.py tests/test_detector.py tests/test_background_routing.py tests/test_folding_fast.py -k 'not physical_isotope_folds_and_limits'
```

The 15 old physical mass-scan cases and complete legacy suite were not rerun;
no full-suite pass is claimed. Initial failed/default-key evidence and successful
logs are retained under `analyses/legend_poster/validation/ionization_acceptance/`.
Existing tests and numerical tolerances were not changed.
