# BG-BENCHMARK-001 implementation report

Implemented and committed on `poster`.

- Commit: `5c62ee81923f75c4a753f8f6e0e8903eb183e7d1`
- Tree: `48d3c6bd4fd55b558a6be0bd4d56f10637e9e794`
- Base: `b2e6a423cb02acc64611d4e52dec047787345b85`

**This is a poster benchmark, not a validated LEGEND-1000 prediction.**

The existing `CompositeBackground` and ER/NR/observed ports now assemble
normalized H-3, tabulated first-forbidden-unique Ar-39, Ge-68/Ga-68, Zn-65 and
Fe-55, pp/B8 solar CEvNS and electron scattering, and an observed residual.
Detector, generic rates/statistics, Ge nuclear datasets, NREFT conventions and
historical notebook remain unchanged. The 26 changed/added files are confined
to backgrounds, source snapshots, focused tests, examples, tools, documentation
and packaging. Per-component detector overrides and existing observed templates
remain available for future measured/MC replacements.

The main corrections are:

- Ar-39 uses a pinned official IAEA BetaShape 2.4 first-forbidden-unique
  calculation, not an allowed beta shape. Its portable table covers 0–565 keV
  and is normalized over full physical support. H-3 retains an explicitly
  approximate allowed spectrum with Coulomb correction and full normalization.
- EC probabilities are isotope- and branch-specific, closed over K/L/M+.
  Reported M and the higher-shell residual remain distinguishable in metadata.
  L and M+ are collapsed to L1/M1 binding energies. Fully contained atomic
  cascades deposit the daughter vacancy energy: Fe-55 K is 6.539 keV rather
  than the approximately 5.90 keV characteristic X-ray energy in the notebook.
- Nuclear excitation and positron branches are counted once per decay, with
  explicit all-escape/full-containment gamma and annihilation outcomes. Rounded
  complete branch sums are explicitly normalized and recorded.
- Closed Ge-68/Ga-68 Bateman evolution and run-averaged activities separate
  exposure/cooldown history from efficiency. Initial daughter activity can be
  zero, secular-equilibrium approximation, or explicitly supplied.
- Physical Ar activity, Ar/Ge mass ratio, transfer probability, veto survival
  and detector acceptance are separate inputs. ER is not NR-quenched; observed
  residuals are not folded again. Discrete lines use detector bin probabilities.

DDEP decay evaluations, ENSDF/NuDat Ar-39 properties, LBNL binding energies,
IAEA BetaShape, original Bahcall spectral tables and B16 solar fluxes are
referenced in the versioned inputs and full `BG-BENCHMARK-001-report.md`.
The starter archive and presentation provided legacy context, not trusted
physical coefficients. All raw table hashes and source-version identities are
retained. The official BetaShape benchmark passed; reference regeneration
reproduced all input files byte for byte. Production requires only NumPy/SciPy;
Matplotlib is an optional plotting dependency.

For the supplied 0.1–15 keV observed window and 365.25 kg-day exposure:

| Scenario | Historical notebook counts | Improved counts |
| --- | ---: | ---: |
| Natural Ge / atmospheric Ar | 1819.58 | 2918.50 |
| Enriched Ge / underground Ar | 153.26 | 110.82 |
| Poster: enriched/UAr, one-year cooldown, 1% veto survival | — | 83.08 |

Ar-39's low-energy shape increases its contribution at the same approximate
full-spectrum benchmark normalization. Cosmogenic counts also change through
run averaging, branching/containment, response acceptance and corrected atomic
energies; these are not isolated single-variable comparisons. The configuration
and per-component counts accompany the figures and CSV.

**Verification:** `python -m pytest -q -W error`: **875 passed in 97.00 s**
(830 existing tests unchanged, 45 new). Tests cover spectral normalization,
K/L/M+ closure, decay branching and Bateman solutions, absolute scattering
checks, narrow pp endpoint integration, routing, exposure scaling, immutable
inputs and integrity failures. The installed wheel contains and loads both the
new background snapshots and existing Ge data. Refining continuum quadrature
changes each integrated poster component by at most 7.313e-7 relative. Exact
patch-tree and bundle restoration checks are recorded separately.

Remaining uncertainties are physical, not resolved by passing tests: notebook
activities/enrichment factors; toy LAr transport with source-shaped Ge deposits;
assumed veto and gamma containment; uncalibrated Gaussian/Fano ER and Lindhard
NR response; shell grouping/atomic escape; approximate minor positron shapes;
pp/B8-only solar coverage, effective-density MSW, atomic ES stepping and
surrogate target masses. Source uncertainties are not propagated. Ar-42/K-42,
Ge-76 double beta decay, full decay simulation, transport and nuisance-aware
likelihood are deliberately absent. Replace components independently as measured
or MC templates become available.
