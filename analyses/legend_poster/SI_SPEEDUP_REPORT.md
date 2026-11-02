# Optional fast SI calculation — implementation report

Branch: **`poster-speedup`**, based on the completed workbench delivery
`f5863ad280e0fc2250b6926d1659920e4f1caf5a`. Later independent changes on
`poster-workbench` are not overwritten or incorporated by this optimization.

The numerical implementation and benchmark source are
**`e76fd698a3f49d24679e4a88bc7662cb212494a9`**, tree
`b8315e0a08bd388c19e3e9190ce543b201d79a3b`. A later input-ownership hardening
copies a supplied breakpoint iterator before checking the explicit SHM match;
it changes neither the numerical kernel nor any benchmark's tuple inputs.
That input hardening is published at
**`f93084236d7770a94e56eed911f47f04142d5e89`**, tree
`d499813d4a577651dc4d390acdbb00fc66551426`.

## Enable and use

In the saved master configuration, set:

```python
working['signal']['si_rate_engine'] = 'factorized_shm'
```

Validate/save in `poster_setup.ipynb`, then run `poster_workbench.ipynb` or the
existing CLI. The default remains `"generic"`; old configurations missing the
field remain generic. No other physical parameter or notebook formula changes.
The active notebooks expose this option.

## Implementation

`halo/shm_eta.py` supplies a vectorized SHM mean inverse speed using the accepted
boosted, hard-truncated distribution. Stable exponential differences handle
zero/tiny boosts and terminal tails. The reference is Freese/Lisanti/Savage,
arXiv:1209.3339v3 Eq. (B4), hard truncation. Narrow integrals use a 16-node
Gaussian evaluation rather than cancellation-prone subtraction.

`si_fast.py` factors the existing O1 rate into the four ordered signed M-response
terms and eta, with the existing nuclear spin and unit factors. It uses the
accepted AFH data, not a Helm model. It requires opt-in and exact zero for every
non-O1 coefficient; a nonzero value, even the smallest representable float,
forces the unchanged generic integrator. Spin/response selection rules and
coefficient-square underflow cannot authorize the fast path.

`SIO1Signal` gains optional `use_fast_si`/`halo` fields and retains direct
physical sigma-to-coefficient normalization and `sigma_upper=mu_upper*sigma_ref`.
Arbitrary signed p/n coefficients, including pure neutron/proton and zero cases,
are supported by the coefficient dispatcher. The finite-ratio sigma_p adapter
retains its existing domain; a pure-neutron model cannot be represented by
sigma_p=0 and a finite cn/cp ratio.

Analysis changes select the optional engine, reuse a mass-specific SI adapter
for limit conversion, and remove exposure/background/threshold-only dependencies
from signal caching where the actual bins/response are unchanged. SI-option
changes do not invalidate non-SI cached spectra. An off-grid threshold still
changes bin integrals; a changed detector response still requires new folds.
All source/data/environment integrity checks remain active.

No generic rate/particle-response formula, frozen convention, AFH dataset,
composition, detector/background/statistical kernel or historical result changed.
The frozen generic standard snapshot remains separate and unchanged.

## Measured performance

These are shared-host measurements, not timing guarantees. Raw scans use the
same prepared models and vectors: 12 masses x 5 isotopes x 24 recoil energies,
two repeats. Construction is excluded equally. Cold folds use the same standard
analysis helper, with separate empty caches, Ge-73 and the configured middle
mass (10 GeV), two repeats; construction/refinement/detector work is included.

| Measurement | Generic median | Factorized median | Speedup |
| --- | --- | --- | --- |
| Raw 1,440-point scan | 231.032 s | 0.109160 s | **2,116.45x** |
| Complete cold single-isotope fold | 35.1908 s | 10.7992 s | **3.2586x** |

Raw samples were 174.354/287.709 s versus 0.121040/0.097280 s. Cold-fold samples
were 38.6101/31.7715 s versus 9.9740/11.6244 s. Contention affects these values;
the raw factor must not be advertised as an end-to-end analysis speedup.

Tracked allocation peaks for a 24-point call were **25,957 bytes generic** and
**7,591 bytes factorized**. This is `tracemalloc` Python/NumPy tracking, excluding
untracked native allocations. Whole-benchmark peak RSS was 105,332 KiB and
includes both paths/imports; it is not a per-algorithm native-RAM comparison.

The separate cold fast calculation of all 35 standard pure-isotope signals took
461.38 s (zero hits, 140 cache misses). Detector folding is now a substantial
remaining cost. This patch deliberately reuses that accepted infrastructure.

## Validation and reproducibility

The numerical benchmark independently compares all 35 O1 pure-isotope folded
spectra and 84 physical SI limits to the immutable generic standard snapshot
(source `a819570f0f0669ab7b286c3dc77fd36a05cf77c4`). Configuration equality is
checked before comparison; changes fail rather than silently comparing unlike
physics. ROI edges agree. Maximum tested-bin relative discrepancy was
`1.9306778398231472e-13`, and maximum relative SI-limit discrepancy was
`4.440892098500626e-16`.

Raw AFH comparisons cover all five isotopes, masses 3/10/100 GeV, equal/unequal/
opposite/global-sign/vanishing p/n couplings, scalar/array/empty inputs and
kinematic boundaries. Every non-SI coefficient/index is checked for exact
fallback, including subnormals; mixed actual generic rates are identical.
Additional checks cover ordered synthetic M terms, 70-digit terminal-tail
oracles, explicit SHM matching, domains/conventions, cache independence and
signed/pure-neutron folded counts and exact 90% CL limits.

The full strict run passed **1099 tests in 4537.31 s**. It was collected before
adding two extra regression functions and the iterator ownership hardening.
After hardening, separate runs passed **175 optimization tests in 154.44 s** and
**219 existing SI/rate/halo tests in 38.45 s**. They overlap: **1101 distinct
cases** are covered; this is not a single 1101-test invocation. Existing tests
were not edited or weakened. `validation/si_speedup/evidence.json` records the
collection/version accounting and exact final log checksums.

The final evidence folder contains benchmark JSON with resolved inputs, actual
source/data hashes, script hashes, reference-manifest hash and both timing
samples, together with regression logs. Reproduce using:

```bash
python -m pytest -q -W error
python tools/benchmark_si_fast.py --output /tmp/wimp-si-benchmark
python tools/benchmark_si_fold.py --output /tmp/wimp-si-fold
```

Detailed interfaces, equations and numerical budgets are in `docs/SI_FAST.md`
and the independent `SI_SPEEDUP_TESTER_HANDOFF.md`. The latter is a handoff,
not a claim of independent Tester approval.

## Remaining assumptions

Only the explicit SHM/O1 factorization is accelerated. Custom PDFs use generic
integration; contradictory fast-SHM declarations fail. Existing calibration,
enrichment, activity, halo and nuclear-source validity uncertainties remain.
No physical range, normalization repair, new coefficient or new detector model
is introduced. The supported mathematical nuclear domain is not a new claim of
physical validity. Prior core audit limitations remain outside this patch.
