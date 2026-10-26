# BG-BENCHMARK-001: low-energy poster benchmark

**Illustrative benchmark, not a validated LEGEND-1000 background prediction.**
This addition replaces notebook background approximations with independently
replaceable ER/NR truth components and an already-observed residual. Detector
ports, generic rate/statistical engines, nuclear-response datasets, frozen NREFT
conventions and the historical notebook are unchanged.

## Reproduce and configure

```sh
python -m pip install -e '.[test,poster]'
python -m pytest -q -W error
python examples/background_benchmark.py --config examples/bg_benchmark_config.json --output bg-results
python tools/background_benchmark/check_convergence.py --output bg-convergence.json
```

The example accepts `--notebook PATH`; by default it reads the repository's
`WIMP_O1_spectrum.ipynb`. It checks its checksum and loads only seven named
function definitions with AST, never executes notebook cells or modifies the
notebook. The last historical function definitions take precedence. A local
compatibility proxy supplies the removed NumPy `trapz` operation. Historical
outputs are already observed and are bin-averaged without further folding.

Outputs: component rates CSV, scenario/configuration/provenance JSON, a 2x2
historical-versus-improved comparison and a poster plot with a threshold-region
zoom, in PNG and vector PDF. CSV rates are events/(kg Ge day keV observed);
figures multiply by 365.25 to display events/(kg Ge Julian year keV observed).

Edit the JSON: initial cosmogenic activities, cooldown and acquisition duration,
Ar activity and Ar/Ge mass ratio, transfer probability, veto survival, nuclear
gamma/annihilation escape assumptions, observed residual, solar inputs, threshold,
resolution, quenching and efficiency are explicit. Initial activities are before
cooldown in decays/(kg Ge day); Ar activity is Bq/(kg Ar). These are not cosmogenic
production rates. `run_duration_day` determines the time-averaged spectrum;
`exposure_kg_day` scales its counts. Changing duration requires rebuilding the
model, whereas scaling exposure with a fixed model scales all counts linearly.

## Public interfaces and data flow

- `AllowedBetaShape(q_keV, daughter_z, daughter_a, beta_minus=True)` is a frozen,
  normalized approximate allowed spectrum, callable for scalar or array kinetic
  energy in keV. Point-Coulomb relativistic Fermi correction; no precision
  exchange/molecular/radiative model. It supplies the finite beta-minus limit at
  zero energy. Production Ar-39 does **not** use this approximation.
- `TabulatedBetaShape(energies_keV, unnormalized_dndE, source_id)` owns immutable
  copies, preserves reported input values and exposes `source_integral`. It
  divides by the full-support piecewise-linear integral, not an analysis-window
  integral, and returns zero outside support. Callers must supply the complete
  physical spectrum; the pinned Ar-39 loader also checks its endpoint.
- `BetaBranch`, `CaptureShell`, `ElectronCaptureBranch` describe branch intensity,
  conditional shell probability and discrete deposited-energy outcomes.
  `ElectronCaptureBranch.lines(..., gamma_escape_probability=...)` returns
  existing `MonoenergeticLine` objects; line folding uses bin probabilities.
- `DecayChain(parent_half_life_day, daughter_half_life_day,
  daughter_production_branch=1)` provides `activities_at` and `mean_activities`.
  An augmented matrix exponential integrates the Bateman equations, including
  equal/near-equal half-lives, without unstable subtractive formulas.
  `mean_activity` supplies the stable single-exponential average.
- `BenchmarkConfig`, `build_benchmark(config, *, er_response, nr_response,
  component_responses=None)` and immutable `BenchmarkModel` assemble components
  on the existing `TruthBackground` / `CompositeBackground` ports.
  Per-component response overrides are supported. `.components` is immutable;
  `.composite` and `.expected_counts(edges, exposure_kg_day)` use existing APIs.
- `backgrounds.neutrinos` exposes `FluxQuadrature`, `SolarSurvivalBenchmark`,
  `SolarScatteringSpectrum` and `helm_squared` for transparent synthetic/physical
  checks. These are separate from the AFH nuclear-response input layer.

Cosmogenic decays and Ar deposition are ER truth: no NR quenching. CEvNS is NR
truth: supplied NR quenching and response. Solar electron scattering is ER truth.
The residual `FlatBackground` is already observed and receives no additional
response or efficiency. Existing `BinnedBackground` can replace any component
with measured/MC observed counts. Activities precede all acceptance factors.

## Physical inputs and approximations

The versioned `backgrounds/data/benchmark_v1/manifest.json` checksums every
input; `benchmark_inputs()` additionally checks a code-pinned manifest hash.
Raw solar tables, BetaShape outputs and source metadata accompany portable
production data. No external decay executable or online retrieval is needed at
runtime. Physical constants and each isotope-specific parameter have references
in `nuclear.json`; no new production target isotope data are introduced.

| Input | Reference and treatment |
| --- | --- |
| H-3 | [DDEP H-3](https://www.lnhb.fr/nuclides/H-3_tables.pdf): 12.312 a, electron endpoint 18.564 keV (distinct from decay Q=18.591 keV); normalized allowed approximation. |
| Ge-68 / Ga-68 | [DDEP Ge-68](https://www.lnhb.fr/nuclides/Ge-68_tables.pdf), [DDEP Ga-68](https://www.lnhb.fr/nuclides/Ga-68_tables.pdf): 270.95 d / 67.83 min, evaluated EC and positron branches with nuclear excitation energies. Exact two-member time evolution after chosen initial daughter activity. |
| Zn-65 | [DDEP Zn-65](https://www.lnhb.fr/nuclides/Zn-65_tables.pdf): 244.01 d; ground/excited EC and 1.421% positron branch, including 1115.549 keV nuclear excitation. |
| Fe-55 | [DDEP Fe-55](https://www.lnhb.fr/nuclides/Fe-55_tables.pdf): 2.75614 a; ground-state EC approximation. Neglect 1.3e-9 excited-state branch and inner bremsstrahlung. |
| EC energy | [LBNL X-Ray Data Booklet Table 1-1](https://xdb.lbl.gov/Section1/Table_1-1a.htm): daughter vacancy energies, not emitted K-alpha photon energies. |
| Ar-39 | ENSDF 39K, J. Chen, NDS 149 (2018), [NuDat decay dataset](https://www.nndc.bnl.gov/nudat3/getdecaydataset.jsp?dsid=39ar+bM+decay+%28268+y%29&nucleus=39K): Q=565(5) keV, half-life 268(8) a, 7/2- to 3/2+, first-forbidden unique. [IAEA BetaShape 2.4](https://github.com/IAEA-NSDDNetwork/BetaShape), X. Mougeot, Appl. Radiat. Isot. 201 (2023) 111018, DOI:10.1016/j.apradiso.2023.111018. |
| Ar activity | AAr 0.964 Bq/kg: [DEAP, EPJC 83, 642 (2023)](https://arxiv.org/abs/2302.14639); UAr 0.00073 Bq/kg: [DarkSide, PRD 93, 081101 (2016)](https://arxiv.org/abs/1510.00702). Source-specific central values, not a LEGEND measurement. |
| Solar flux | B16-GS98 totals pp=5.98e10 and B8=5.46e6 cm^-2 s^-1, [Vinyoles et al., ApJ 835, 202 (2017), Table 6](https://arxiv.org/abs/1611.09867). |
| Solar shape | [Bahcall pp table](https://www.sns.ias.edu/~jnb/SNdata/Export/PPenergyspectrum/ppenergytab), PRC 56, 3391 (1997), Table IV; [Bahcall et al. B8 central table](https://www.sns.ias.edu/~jnb/SNdata/Export/B8spectrum/b8spectrum.txt), PRC 54, 411 (1996). Original files retained; explicitly normalize their rounded trapezoidal integrals 1.000003343 and 0.999986544. |
| Weak/scattering | GF=1.1663787e-5 GeV^-2; sin^2 thetaW=0.23863 low-energy effective benchmark, [PDG 2023 electroweak review](https://pdgweb.lbl.gov/2023/reviews/rpp2023-rev-standard-model.pdf). Leading Freedman CEvNS, Helm form factor (Lewin & Smith, Astropart. Phys. 6, 87 (1996)); tree-level neutrino-electron scattering with explicitly supplied constants. |

EC conditional probabilities are isotope- and branch-specific. DDEP K, L and
reported M probabilities are retained separately. For a closed benchmark shell
model, **M+** groups the reported M and the residual higher-shell probability;
raw M and residual are recorded rather than mislabelling M+ as pure M. K, L and
M+ sum to unity. Representative ground-state (K,L,M+) values are Ge-68
(0.8639,0.1150,0.0211), Fe-55 (0.88044,0.10238,0.01718), and Zn-65
(0.8853,0.0977,0.0170). Other EC branches use their own evaluated probabilities.
L and M+ are collapsed to L1 and M1 vacancy energies. This neglects sub-shell
splitting, solid-state shifts and atomic-cascade escape.

Fully contained atomic cascades deposit daughter binding energy: Fe-55's K
feature is **6.539 keV**, not the notebook's approximately 5.90 keV Mn K-alpha
photon. Ge-68 and Zn-65 K energies are 10.367 and 8.979 keV. Each nuclear gamma
cascade is approximated as wholly escaping or wholly contained; positron
annihilation likewise deposits either zero or 1022 keV, independently. Beta
kinetic energy remains contained. No partial Compton energy transport, IC atomic
cascades, angular correlations or delayed Ga coincidence tagging is simulated.
Escape probabilities are adjustable geometry assumptions, not branching ratios.

Rounded complete branch totals are 1.00001 for Zn-65 and 0.9999416 for Ga-68.
These reported values and explicit normalization factors are retained in
provenance, then branches are divided by those totals. We use a consistent DDEP
set; newer ENSDF evaluations with different probabilities/denominators are
not mixed into it. Minor nonallowed Ga/Zn branches are treated with an allowed
positron shape; small fractions and large containment offsets limit their
impact in this window. This is an approximation, not precision decay modelling.

Solar coverage is **pp and B8 only**, not the full solar spectrum. Be-7, pep,
CNO and hep are absent. ES uses an energy-dependent adiabatic/decohered MSW
survival probability at one configurable effective production density; no
production-region averaging or Earth matter treatment. Atomic ES uses a Ge
shell stepping approximation, four valence electrons treated as free. CEvNS
uses the leading heavy-target endpoint and all active flavors. The effective
A=72.63/76 targets, masses A*0.93149410242 GeV and target counts NA*1000/A are
explicit **surrogates**, not precise isotope masses or enriched compositions.
This approximation remains isolated from the production nuclear-response data.
The assumed Gaussian/Fano response and Lindhard extrapolation to few-eV pp
recoils are not calibrated. Mathematical full support is implemented; no
experimental validity range is claimed.

## Benchmark assumptions and comparison with the notebook

The starter archive was used as a source inventory and independent beta-shape
check, not promoted to authoritative data. The presentation's approximate
component comparison is retained. H-3 77/year, natural-Ge Ge-68 600/year,
Zn-65 60/year and Fe-55 12/year, and enrichment divisors 14, 2.5 and 1.33 are
**notebook normalization assumptions**, not accepted production measurements.
Cosmogenic activation literature establishes their strong exposure-history
dependence (e.g. GERDA detector production, EPJC 74, 3253 (2014),
[arXiv:1410.0853](https://arxiv.org/abs/1410.0853);
Wei, Mei & Zhang, Astropart. Phys. 96 (2017),
[arXiv:1706.05324](https://arxiv.org/abs/1706.05324)). No production rate is
silently substituted for a decay activity.

The LAr benchmark assumes 20 kg Ar/kg Ge and transfer probability
0.00014230663956447948, calibrated to the notebook's full-spectrum Ar
normalization (86583.70384812194 events/kg-year before response for its AAr
reference). The conditional energy deposited in Ge is assumed to have the
source beta shape. **This is a normalization/transport surrogate**, not an
LAr electron-transport prediction. Transfer and veto survival multiply the
physical Ar decay activity independently; detector efficiency is applied once
by the assigned ER response. AAr/UAr depletion now follows the two cited
activity measurements (ratio approximately 1321), not the notebook factor 1400.

Comparison profiles have zero cooldown, one-year acquisition and veto survival
1. The poster profile uses enriched/UAr settings, one year cooldown, one year
acquisition and **assumed 1% veto survival**. Default nuclear gamma and
annihilation escape probabilities are zero (complete containment). Initial
Ga-68 activity equals initial Ge-68 activity as a secular-equilibrium
approximation; the exact Bateman solution then evolves it. Zero or explicit
initial Ga activity are supported. No independently counted equilibrium daughter
is added on top of a parent line.

All profiles use a 0.1 keV analysis threshold, 0.05 keV noise sigma plus Fano
term, constant 0.9 efficiency and configurable Lindhard NR quenching. Residual
rates 0.1/0.115 events/(kg-year keV observed) already include their acceptance.
For 365.25 kg-day exposure in the **0.1–15 keV observed** window:

| Scenario | Historical notebook | Improved benchmark |
| --- | ---: | ---: |
| natural Ge / AAr | 1819.58 | 2918.50 |
| enriched Ge / UAr | 153.26 | 110.82 |
| enriched/UAr poster, cooldown + veto assumptions | not a historical profile | 83.08 |

These changes combine physics and explicit acquisition/response choices, not a
controlled single-variable physical correction. Ar-39 changes most: its
forbidden/Coulomb shape increases low-energy counts (AAr 1081.33 to 2443.72 in
this window at the chosen normalization). EC totals decrease through decay
averaging, branching/containment and acceptance; Ga-68 becomes explicit; M+
appears below about 0.2 keV; Fe-55 K shifts upward. H-3 is normalized over its
full endpoint, then decayed/folded, rather than normalized over the plotting
window. Solar spectra now have independently sourced flux/shape inputs; their
coverage and detector limitations must accompany any poster use.

## Regeneration and numerical verification

```sh
# Optional Linux x86-64 reference-tool regeneration; runtime is portable.
# Python >=3.12 needed for the regeneration tar extraction filter.
python tools/background_benchmark/regenerate.py --work /tmp/bg-sources --output /tmp/bg-inputs
# Compare /tmp/bg-inputs with src/wimp_legend/backgrounds/data/benchmark_v1.
```

The script pins BetaShape repository commit
`22b8691dd3ea34e6b69fcb93167ab872673ffd5b`, archive SHA-256
`fa0a90094084b417c3852bcc93adcadba097d6d05187516c9aae1241595ad30e`,
and original solar-file hashes. It verifies extracted package members, runs the
official supplied benchmark and checks for reference differences, then runs
`betashape Ar39.txt -csv myEstep=0.1`. Its minimal ENSDF input is an explicitly
labelled transcription of evaluated Ar-39 quantum numbers, endpoint and branch,
not an original evaluator download. Raw output includes numerical screening,
exchange, overlap, radiative and first-forbidden-unique corrections. Source
point uncertainties are retained but not propagated; requested 0.1 keV output
interpolates the tool's default calculation grid. All ten input files and the
manifest regenerate byte for byte in the verified reference environment.

Original delivery validation: **875 passed**, `python -m pytest -q -W error` (830 existing tests
unchanged, 45 new). New checks cover full beta normalization, independently
known spectral limits, tabulated input ownership, K/L/M+ sums, branch/escape
conservation, Bateman solutions/time averages, solar flux and scattering
normalization, narrow pp endpoints, ER/NR/observed routing, exposure scaling,
response overrides, observed-template replacement and checksum tampering.

The plotting-only fold uses existing detector bin probabilities and exact
point-line probabilities. Continuum fixed Gauss-Legendre quadrature uses
source breakpoints and bin crossings. Additional pp NR breakpoints prevent
ordinary adaptive integration from missing its very narrow few-eV support.
The example provides a rigorous Gaussian sigma bound for an explicit 12-sigma
continuum truncation (missing probability <= 2*Phi(-12) per event); generic
responses have no such truncation. GL6/0.05 keV compared to GL10/0.025 keV
changes each component's integrated poster-window count by at most
**7.32e-7 relative**. Detailed per-bin errors and refined counts are in the
convergence JSON. This establishes numerical convergence of the chosen model,
not independent physical validation.

## Default public-count integration audit follow-up

The original 875-test result did not cover the default-policy Ar39 and solar
ES public folds. CODE_AUDIT_007 reproduced exceptions in those paths while
independently accepting the separate poster plotting calculation. The factory
now supplies all Ar39 interpolation knots and all fixed-flux ES/CEvNS recoil
endpoints, and the shared integration policy batches dense breakpoint sets
without loosening tolerances, dropping source support or suppressing warnings.
The original physical inputs and approximations remain unchanged. New
regressions exercise the ordinary `expected_counts` methods with the default
policy, including full-support and Gaussian folds against analytic/CDF
references. See [CODE_AUDIT_007_RESPONSE.md](CODE_AUDIT_007_RESPONSE.md) for the
correction and its independent-review boundary.

## Replacement priorities and independent review

1. Replace initial activities with measured exposure/cooldown histories and
   actual detector compositions; audit DDEP/ENSDF evaluation choices and minor
   branch/shell approximations against the desired evaluation release.
2. Replace LAr transfer and gamma/annihilation escape models with measured or MC
   deposited-energy templates; establish veto/acceptance correlations.
3. Supply calibrated low-energy ER/NR response, atomic escape, thresholds and
   quenching; treat source uncertainties outside this fixed benchmark.
4. Add complete solar flux coverage and production-region survival/atomic
   ionization calculations before calling the neutrino floor quantitative.

Optional Ar-42/K-42 and Ge-76 2nu-beta-beta are **not implemented**. This task
does not supply decay chains beyond Ge-68/Ga-68, a transport simulation, an
empirical calibration or a nuisance-aware likelihood. No change to frozen
NREFT conventions is needed or made.
