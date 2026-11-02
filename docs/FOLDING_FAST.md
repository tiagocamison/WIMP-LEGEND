# Optional vectorized detector folding

Set these independent options in the master `analyses/legend_poster/config.json`
(or edit and save them from `notebooks/poster_setup.ipynb`):

```python
working['signal']['si_rate_engine'] = 'factorized_shm'
working['numerical']['folding_backend'] = 'vectorized'
```

Run `poster_workbench.ipynb` normally. Both defaults remain the existing
reference implementations (`generic` SI rates, `reference` folding). An old
configuration without `folding_backend` uses `reference`. The folding option
also works with generic NREFT spectra; it does not select or alter a rate kernel.

## Scope and mathematics

The workbench still delegates continuum integration to
`tools/background_benchmark/folding.py::quadrature_counts`. Gauss-Legendre
orders, truth steps, bin crossings, source/interpolation breakpoints, support
clipping, Gaussian tail bound and coarse/fine convergence tests are unchanged.
It still integrates *bin probabilities*, not sampled observed-energy densities.

`VectorizedGaussianBins` computes the same Gaussian probability

\[
p_j(E)=\epsilon_{\rm obs}\left[\Phi(z_{j+1})-\Phi(z_j)\right],\quad
z_j=(E_j-\bar E(E))/\sigma(\bar E(E)).
\]

For positive lower endpoints it instead uses
`Phi(-z_j) - Phi(-z_{j+1})`, exactly as the reference, to avoid cancellation.
NR uses the existing quenching model; ER mean energy is the ER truth energy.
The existing truth efficiency multiplies the spectrum in quadrature, outside
this kernel. There is no finite-window renormalization, extra efficiency,
interpolated observed count or change in physical units.

Only concrete built-in ER/NR responses with the built-in Gaussian resolution
and absent/constant observed efficiency use the array calculation. Custom
resolution/response ports retain their own `expected_counts` implementation;
nonconstant observed efficiencies use the original probability integration.
Point lines always delegate the original, unpruned `expected_line_counts` API.
Already-observed spectra/counts continue to bypass detector processing.

The existing continuum pruning removes bins outside the same batch-union
12-sigma window (or the configured tail bound). It is **not** a new row-wise
truncation. Temporary CDF arrays use at most 128 truth rows at a time. There is
no JIT, fastmath, thread/process pool or new dependency; compilation overhead
is zero. The existing configured BLAS thread limit still applies.

## Cache and reuse

The same content-keyed numerical cache is retained, including source/data and
environment fingerprints, checksum validation and owned finite arrays.
`folding_backend` is included in both background and folded-signal cache keys.
Switching backends cannot read a reference fold as an optimized fold (or vice
versa). Truth-rate keys exclude the folding backend, detector, binning,
exposure and background. Truth-rate tables can therefore be reused when only
folding changes. The existing per-run isotope spectra are reused across
composition/scenario/exposure combinations where their response and bins match.

As before, exposure and styling changes reuse unit-exposure numerical counts;
an on-grid threshold can reuse folds, while a threshold that inserts a new
edge correctly requires bin reintegration. Changing response noise, efficiency,
quenching or actual bin edges invalidates the relevant folds. No approximate
kernel matching, global mutable response cache or unrelated-background key is
introduced for truth signals.

## Verification and timing

See `analyses/legend_poster/FOLDING_SPEEDUP_REPORT.md` and
`FOLDING_SPEEDUP_TESTER_HANDOFF.md` for exact tested source identity and evidence.
Reproduce the complete matched cold/warm dense O1 scan with:

```bash
python tools/benchmark_folding.py --output /tmp/fold-scan --masses 24
```

Use a new output directory each time. Each backend gets a separate process and
empty cache before its cold run, then the same configuration/cache for its warm
run. Both use the fast SI engine, standard integration settings, all configured
scenarios and exposures, and the historical comparison. The 24-point mass grid
is log-spaced from 3 to 100 GeV; all five natural Ge isotopes contribute.

`--masses 0` retains the seven standard masses. For a complete-analysis profile,
add `--profile-analysis` and use a separate output directory. Profiling overhead
is included in those timings; do not mix them with the unprofiled benchmark.
The tool checks all saved signals/backgrounds, ROI boundaries, limits and
convergence. Peak RSS is measured per backend subprocess, including imports and
both cold/warm calls; it is not a claim about a native-kernel-only allocation.

The optimization makes no new claim about LEGEND calibration, nuclear-data
validity, or independent scientific acceptance.
