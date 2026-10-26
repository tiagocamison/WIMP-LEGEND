# LEGEND poster analysis workbench

This is a **provisional benchmark**, not a validated LEGEND-1000 prediction.
The working input is [`config.json`](config.json). The two notebooks and CLI
call the same Python orchestration; neither contains a second set of physical
parameters. The accepted package physics, nuclear datasets and conventions are
unchanged.

## Start here

From the repository root, in a fresh Python 3.12 environment:

```bash
python -m venv .venv
.venv/bin/python -m pip install -r analyses/legend_poster/requirements.txt
.venv/bin/python -m pip install -e '.[test,workbench]'
.venv/bin/python -m ipykernel install --user --name wimp-poster
```

Select that kernel for both notebooks. A Jupyter frontend must be installed
separately if your editor does not provide one. Python dependencies for the
tested environment are pinned in `requirements.txt`; the package's general
dependency bounds are not tightened.

1. Open `notebooks/poster_setup.ipynb`. Inspect the target, detector, signal,
   background, analysis and plotting dictionaries. Edit `working`, validate it,
   and run the explicit `save_config` cell. Examples are commented out so merely
   executing the notebook does not change physical assumptions.
2. Open `notebooks/poster_workbench.ipynb`. It reads the **saved** configuration,
   computes or retrieves numerical outputs, plots them, and displays limits,
   counts, ROI endpoints and convergence diagnostics.
3. Start with the smoke profile to check installation. It evaluates all three
   operators and all five natural-Ge isotopes but only two masses, with explicitly
   coarse numerical tolerances. Do not use smoke figures on the poster.
4. Use the standard profile for the configured seven-mass grid. The first run
   is more expensive; subsequent exposure/style changes reuse the rates and
   folded spectra. Refine the mass grid for a smooth final exclusion curve.

The environment variables `WIMP_POSTER_CONFIG`, `WIMP_POSTER_PROFILE`,
`WIMP_POSTER_CACHE`, and `WIMP_POSTER_RUNS` select a saved configuration, profile,
cache and output directory. These select inputs/locations, not hidden physics.
The default profile is specified in the master file.

Equivalent commands, from the repository root:

```bash
python -m analyses.legend_poster --profile smoke
python -m analyses.legend_poster --profile standard
python -m analyses.legend_poster --plot-only analyses/legend_poster/runs/RUN_ID
python -m pytest -q -W error
```

`--plot-only` reads the saved numerical run and the current master plotting
section; it does not run a detector fold or rate calculation. All figures are
exported as vector PDF and PNG. Run/cache directories are deliberately ignored
by Git; they can be regenerated.

## Configuration and scientific conventions

`load_config`, `validate_config`, `save_config`, and `resolve_profile` are public
through `analyses.legend_poster`. Loading/validation own their dictionaries;
save validates before atomically replacing the file. Invalid configurations
fail rather than silently renormalizing fractions or converting conventions.

* **Targets:** all five Ge isotope records are explicit. Natural number
  abundances and atomic masses are a NIST SRD 144 snapshot with stated
  uncertainties and access date, not a claim to the latest CIAAW evaluation.
  The 92% Ge-76 / 8% Ge-74 number-fraction preset is **hypothetical**, not an
  adopted LEGEND enrichment. Replace it with a sourced measurement when available.
* **Masses:** kinematics uses bare nuclear masses derived from sourced atomic
  masses, CODATA 2022 conversion/electron mass and the AME2020 approximate total
  electron-binding fit. That small correction is labeled approximate; no mass
  uncertainty propagation is performed. `A * 0.9315` is not used as isotope data.
* **Fractions and kilograms:** `TargetComposition` mass fractions retain their
  existing *bare nuclear ensemble* semantics. `isotope_weights` explicitly
  converts to number fractions, then computes neutral-atomic bulk mass fractions
  and nuclei per kg. A pure-isotope rate uses nuclei/kg of that neutral isotope;
  its mixture weight is the neutral-atomic bulk mass fraction. Thus
  `bulk_fraction_i * pure_nuclei_per_kg_i == mixture_nuclei_per_kg_i`.
  Electrons contribute detector mass, not recoil kinematics. Number and mass
  input presets denoting the same ensemble give the same final rates.
* **Detector:** NR uses the existing Lindhard quenching plus Gaussian response;
  ER uses the existing ER response without quenching. Constant observed-energy
  efficiency is applied once to truth components. Threshold is an analysis cut
  on observed bins, not an extra efficiency. Noise, Fano and efficiency values
  are uncalibrated benchmark assumptions. Components may select independent
  response presets with `components.NAME.response`; observed residuals cannot.
* **Backgrounds:** the complete existing `BenchmarkModel` is built, including
  H3, Ge68/Ga68, Zn65, Fe55, Ar39, solar CEvNS, solar ES and the observed flat
  residual. Its source-data files/checksums and physical kernels are reused.
  Per-component `enabled`, `scale`, and response overrides are explicit. Activities
  are physical decays/(kg day), distinct from LAr transfer/veto or containment.
  Cooldown and acquisition duration determine the averaged physical activities;
  an exposure scan scales those *same averaged rates*. Change acquisition/cooldown
  explicitly to model a different time history. No new decay/transport model is
  introduced. The accepted benchmark has pp/B8 solar inputs, not a complete
  Be7/pep/CNO/hep inventory, and no Ar42 or Ge-76 double-beta component.
* **Solar mixture:** CEvNS delegates to the existing isotope-specific kernels
  with the configured nuclear masses, A and nuclei/kg. ES uses the exact total
  atom count and common Ge Z/shell kernel; it does not need an effective isotope.
  `effective_isotope` remains an explicitly labeled alternative. Effective
  factory bookkeeping in the inherited background provenance is superseded by
  `workbench_solar_target` in mixture mode.
* **Signal:** accepted AFH germanium datasets, oscillator prescription and
  sub-per-mille source closure deficits are unchanged. Halo velocities are km/s
  at the explicit SHM constructor boundary. The mass grid, density, spin, signs,
  neutron/proton ratio and numerical nucleon reference mass are explicit.
  Frozen `anand2014_prc89_065501_half_isospin_v1` is required exactly.
  The same `fn_over_fp` means `c_n/c_p` for each selected single operator.
* **O1:** the existing physical `SIO1Signal` adapter defines the proton SI
  cross section, including the existing GeV^-2-to-cm^2 constant. The returned
  cross-section upper limit is `mu_upper * sigma_p_ref_cm2`.
* **O8/O11:** one operator at a time, dimensionful proton/neutron coefficients
  in GeV^-2; `c_p = proton_sign * coefficient_ref_GeV_minus2` and
  `c_n = fn_over_fp * c_p`. Plots show `|c_p|_upper = sqrt(mu_upper) * |c_p_ref|`,
  **not** SI cross sections. Default spin is 1/2. No automatic convention/basis
  adapter is present. `m_N_reference_GeV` is explicit configuration/provenance;
  the analytic particle kernel receives Q through the existing rate assembly.
* **Statistics:** the default is the existing exact Poisson-CLs counting
  construction with its documented expected-background pseudo-observation.
  The threshold-anchored ROI ends at the first bin containing at least 95% of
  the signal **inside the configured observed analysis window**. It is not a
  promise to contain 95% of all unobserved recoils. Counts, actual retained
  fraction, pseudo-observation, CL and construction are saved. Optional
  `binned_asimov` uses the full above-threshold window and is separately labeled;
  it is not the historical method or a nuisance-aware likelihood.

The historical `WIMP_O1_spectrum.ipynb` background is an already-observed
reference. It is integrated anew in the requested observed bins without another
response or efficiency. The paired revised/historical limits use the **same new
signal, exposure, ROI and statistics**. Historical background response/selection
remain frozen: this is not a recalibration of the old background under arbitrary
new detector assumptions. The pilot notebook is preserved byte-for-byte at
`notebooks/historical/poster_workbench_pilot.ipynb`; original assets remain in
`notebooks/poster_workbench_assets/` and never feed production sensitivities.

## Numerics, saved results and caching

`run_analysis(config, profile=..., cache_dir=..., output_dir=...)` returns a
checksum-validated `load_run` result. Each run writes:

* `backgrounds.npz`: integrated events/(kg day) in every nonuniform observed bin,
  each component and total, plus the historical comparison when selected.
* `signals.npz`: reference integrated events/(kg day), both pure-isotope rates
  and combined scenario rates; `limits.csv`: mass-dependent limits and counts.
* `convergence.json`: component coarse/fine quadrature checks and isotope/signal
  truth-grid refinement checks; failed convergence aborts rather than certifies.
* `manifest.json`: full resolved master configuration, Git commit/dirty state,
  actual core/analysis/helper/data SHA-256 hashes, environment versions, derived
  masses/weights, statistics, numerical settings and output checksums.
* `figures/plot_manifest.json`: numerical-manifest hash, current plotting settings,
  external-data provenance and figure hashes. Plots reload saved numerical data.

Near threshold the default bins are 0.01 keV; higher bins are 0.05 keV. All
observed bins are integrated afresh against detector probabilities. Coarse
observed counts are **never interpolated**. PCHIP is used only for nonnegative
truth recoil rates evaluated by the accepted core; truth grids are adaptively
doubled and checked. Negative core rates raise rather than hide a source issue.
Quadrature order/step are independently refined. Convergence is measured for
total counts and bins above a configurable fraction of the largest bin.

For continuum speed, the existing quadrature helper and a local Gaussian-bin
wrapper delegate to the core CDF only in the relevant bins. Source/bin pruning
uses an explicit sigma distance (default 12); the combined omitted probability
is bounded by `erfc(tail_sigma/sqrt(2))` per event and validated against the
configured maximum. Discrete lines always use the exact unpruned core line-bin
probabilities. The approximation is numerical, not a new detector calibration.

Cache keys include resolved relevant inputs, actual implementation/data hashes
and Python/NumPy/SciPy versions. Only finite numerical arrays are stored, without
pickle; their checksums are checked on load. Raw truth rates can survive detector
changes; folded counts can survive exposure/style changes. Failed/partial runs
can reuse completed kernels. Plotting code is excluded from physics cache keys.
Do not delete the cache merely to change a legend or figure range.

## External exclusions

`data/si_limits/manifest.json` freezes individual observed limits and projections
from cajohare/DirectDetectionPlots, with upstream commit, paths, original headers,
raw SHA-256/Git blob hashes, units, confidence levels and references. Observed and
projection curves have different line styles. Their astrophysical/statistical
assumptions are not transformed to ours. These snapshots may be outdated in
2026. CDMSLite has multiple references in its upstream header and must not be
relabeled as a uniquely reconstructed release. The external notebook and plotting
utilities were inspected for organization; their code is not imported or executed.

See `IMPLEMENTATION_REPORT.md` and `TESTER_HANDOFF.md` for validation and remaining
physical limitations. No core physics, historical notebook, nuclear dataset or
detector/statistics module is modified by this workbench.
