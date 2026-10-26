# Independent-tester handoff: BG-BENCHMARK-001

Review target: commit `5c62ee81923f75c4a753f8f6e0e8903eb183e7d1`, tree
`48d3c6bd4fd55b558a6be0bd4d56f10637e9e794`, branch `poster`.
This handoff specifies review work; it does not claim an independent reviewer
has certified the benchmark.

## Restore and run

```sh
git clone --branch poster WIMP-LEGEND-poster.bundle WIMP-LEGEND
cd WIMP-LEGEND
python -m pip install -e '.[test,poster]'
python -m pytest -q -W error
python examples/background_benchmark.py --config examples/bg_benchmark_config.json --output bg-results
python tools/background_benchmark/check_convergence.py --output bg-convergence.json
```

Alternatively unpack `WIMP-LEGEND-source.zip`. The unified diff applies to
base `b2e6a423cb02acc64611d4e52dec047787345b85`. Archival fixed-column ENSDF
and original source tables intentionally preserve trailing whitespace. Use
`git apply --whitespace=nowarn BG-BENCHMARK-001.diff`; do not auto-fix source
bytes. Expected complete suite: **875 tests**. See `identity.json` for versions,
protected hashes, and `validation/` for actual run logs and numerical checks.

## Public interfaces to probe

`wimp_legend.backgrounds` exports `AllowedBetaShape`, `TabulatedBetaShape`,
`BetaBranch`, `CaptureShell`, `ElectronCaptureBranch`, `DecayChain`,
`mean_activity`, `BenchmarkConfig`, `BenchmarkModel`, `build_benchmark` plus
unchanged existing background ports. `backgrounds.neutrinos` exposes
`FluxQuadrature`, `SolarSurvivalBenchmark`, `SolarScatteringSpectrum`,
`helm_squared`. `backgrounds.benchmark.benchmark_inputs()` returns fresh
independent metadata and an immutable tabulated Ar-39 shape after hash checks.

| Property | Independent expectation |
| --- | --- |
| Beta PDFs | Integral over complete physical support is one; nonnegative; scalar/vector agreement; zero outside support; changing analysis bins does not renormalize. |
| Ar-39 physics | 7/2- to 3/2+ is 1U; endpoint 565 keV; official tool mean consistent with evaluated approximately 218.8(2.1) keV; no allowed-shape fallback. |
| EC | Conditional K/L/M+ probabilities sum one per branch. Sum of all gamma/annihilation outcomes equals that branch's normalized probability for every escape parameter. Daughter vacancy rather than K-alpha energy is deposited. |
| Chains | dAp/dt=-lambdaP Ap; dAd/dt=b lambdaD Ap-lambdaD Ad. Recover analytic distinct-half-life solution, equal-half-life limit, zero daughter growth and secular asymptotics. Run means are integrals/duration, not endpoint activities. |
| Flux | Sum quadrature weights equals configured total flux cm^-2 s^-1. Pee in [0,1], correct vacuum/high-energy limits. CEvNS all active flavors, ES separately mixes flavors. |
| Scattering | Free-electron integrated ES agrees with analytic integral; leading CEvNS total at F=1 agrees with GF² Qw² Eν²/(4pi). Convert GeV^-2 once and dT GeV to keV once. |
| Routing | ER bypasses NR quenching; CEvNS receives assigned NR response. Observed backgrounds/counts never receive another efficiency. Lines use exact bin probabilities, not quadrature of delta distributions. |
| Scaling | For a fixed model, counts are linear in exposure. Ar source activity, mass ratio, transfer and veto multiply once each. Changing acquisition duration rebuilds averaged activity; it is distinct from scaling exposure. |
| Replacement | Per-component response override affects only that component; substituting BinnedBackground requires no physics-engine edits. |
| Ownership | Caller array/dict edits do not alter frozen stored inputs; returned metadata edits do not mutate later loads; manifest/table corruption is rejected. |

## Reproduce physical sources

Optional Linux x86-64 / Python >=3.12:

```sh
python tools/background_benchmark/regenerate.py --work /tmp/bg-sources --output /tmp/bg-inputs
diff -rq /tmp/bg-inputs src/wimp_legend/backgrounds/data/benchmark_v1
```

No BetaShape executable is needed for ordinary package use. The regeneration
script checks official package/source hashes, supplied benchmark reference
results and portable exports. Independently audit the minimal ENSDF
transcription, default corrections and branch data; a passing regeneration
only proves consistency with this recorded tool run.

Review DDEP evaluation versions and PM-to-M+ treatment; raw rounded totals
and normalization factors must remain visible. Solar tables use full-source
normalization, four-point quadrature on each original spectral interval; pp NR
support is narrow and must not disappear in adaptive folding. Compare both
integrated and per-bin values in `validation/convergence.json`.

## Physical caveats to carry into poster use

The LAr transfer probability is a normalization/transport assumption with a
source-beta-shaped conditional Ge deposit, not an electron transport simulation.
The poster's 1% veto survival is assumed. Atomic cascades are fully contained;
gamma/annihilation only have all-or-nothing outcomes. No measured enrichment
composition or activity history is supplied. Solar coverage is pp/B8 only;
other fluxes, production-region survival averaging and accurate atomic response
are absent. Effective target masses and low-energy Lindhard/Gaussian response
are unvalidated approximations. Independent source/physics checks are needed
before drawing quantitative LEGEND sensitivity conclusions.

The full normative NREFT document, production Ge response datasets and historic
notebook must remain byte-identical to the base. No automatic convention or
nuclear-response transformation is part of this task.
