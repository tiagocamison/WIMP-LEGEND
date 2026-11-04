# Technical tester handoff

Scientific development candidate; independent acceptance is external.

Apply supplied patch to feature/pre-smearing-acceptance and validate affected interfaces

## Source and authorization

```json
{
  "source": {
    "commit": "5e744a2dd3bff513abbe82cbbadee935768ee0c3",
    "tree": "f380defcbcb09edb158b893505250f44a89fb1af",
    "dirty": false,
    "inputs": {}
  },
  "authority": {
    "source": "User 2026-11-04: apply ionization_acceptance.patch, check branch and push upon completion",
    "scope": "Patch application, compatibility fixes, focused regression testing and branch push",
    "publish": true,
    "merge": false
  }
}
```

## Status

Implementation: complete; internal verification: passed; independent acceptance: pending.

## Assumptions and criteria

- No calibration data added; acceptance disabled by default
- No configured environment_status tool or network-policy snapshot available; public Git via ordinary configured proxy
- C1 [pass]: Acceptance patch and affected detector/background regressions pass — 121 affected tests passed with warnings as errors, including all 19 ionization tests. Existing detector/routing/folding regressions passed after restoring the branch missing explicit reference default. Fifteen old physical mass-scan cases intentionally not rerun. (runs: run-0002)

## Execution evidence

- run-0001 [failed], exit 1, charged 54.534 s: `/workspace/scratch/ebd48d783ed1/acceptance-evidence/run-0001/run.json`
  argv: ["../acceptance-env/bin/python", "-m", "pytest", "-q", "-W", "error", "tests/test_ionization_acceptance.py", "tests/test_detector.py", "tests/test_background_routing.py", "tests/test_folding_fast.py", "-k", "not physical_isotope_folds_and_limits"]
  persistence: {"location": "https://github.com/tiagocamison/WIMP-LEGEND/tree/53f7b4ae47e14ec924293fc5f7ff9ef88470564e/analyses/legend_poster/validation/ionization_acceptance/failed-initial/", "manifest_sha256": "25c493e23d0962139c9a624e5fc47b0261d8560cf09b776556d898f1d3193b6c"}
- run-0002 [passed], exit 0, charged 63.846 s: `/workspace/scratch/ebd48d783ed1/acceptance-evidence/run-0002/run.json`
  argv: ["../acceptance-env/bin/python", "-m", "pytest", "-q", "-W", "error", "tests/test_ionization_acceptance.py", "tests/test_detector.py", "tests/test_background_routing.py", "tests/test_folding_fast.py", "-k", "not physical_isotope_folds_and_limits"]
  persistence: {"location": "local only"}

Budget: {"authorized": true, "provenance": "User authorized check/apply/push; focused tests bounded at 20 minutes, no dense scans or full legacy suite", "wall_seconds": 1200, "max_attempts": 8, "ram_gib": 3, "max_concurrent_commands": 1}
Consumed wall seconds: 118.380

## Limitations and blockers

- None recorded.

## Next action

User/Tester: Pull feature/pre-smearing-acceptance and choose none/hard/erf in the saved detector configuration.
