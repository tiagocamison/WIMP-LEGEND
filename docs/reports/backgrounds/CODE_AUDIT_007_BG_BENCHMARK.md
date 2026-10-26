# BG-BENCHMARK-001 — independent scientific audit

Audit date: 2026-10-25. **Poster plotting path: PASS WITH CONDITIONS. Default public count API: FAIL.** The delivered plotting model is independently reproduced, but normal calls to `expected_counts` can fail for Ar39 and solar electron scattering with the default integration policy. The blanket claim that only physical uncertainties remain is therefore too strong.

No production code, repository tests or previous audit assertions were changed. The strict failing tests are delivered as failures.

| Area | Verdict | Scope |
| --- | --- | --- |
| Source identity and snapshot integrity | **PASS** | Exact commit/tree restored from delivery bundle; diff independently reproduces tree; protected source unchanged. |
| Physical input provenance | **PASS WITH CONDITIONS** | Ar39 calculation independently rerun; dominant evaluated decay/atomic/flux inputs checked. Complete primary PDF transcription audit remains incomplete. |
| Beta spectra | **PASS WITH CONDITIONS** | Normalization, support and point-Coulomb approximation agree with high-precision references; official Ar39 first-forbidden-unique calculation reproduced. H3 and minor positron spectra remain approximations. |
| EC branches and decay chains | **PASS WITH CONDITIONS** | Correct daughter vacancy energies, conditional shell bookkeeping, cascade probability conservation and Bateman means; unresolved subshell/cascade/transport approximations remain. |
| Solar scattering | **PASS WITH CONDITIONS** | Absolute cross sections, units, MSW limits, flux normalization and threshold integrals pass. pp/B8 coverage and detector/target surrogates do not constitute a neutrino-floor calculation. |
| ER/NR/observed routing and scaling | **PASS** | Independent line probabilities and poster integrals confirm mapping and acceptance; overrides, exposure and observed-template replacement tested. |
| Default public count integration | **FAIL** | Reproducible numerical failures for ordinary poster-window Ar39 and ES counts; see BG-AUDIT-001/002. |
| Supplied plotting example | **PASS WITH CONDITIONS** | Works through its separate numerical folding helper; independent component counts agree within approximately 8e-7 relative. Conditions are benchmark physics assumptions, not evidence of LEGEND sensitivity. |

## Exact tested source

| Item | Identity |
| --- | --- |
| Commit | **`5c62ee81923f75c4a753f8f6e0e8903eb183e7d1`** |
| Tree | **`48d3c6bd4fd55b558a6be0bd4d56f10637e9e794`** |
| Base commit / tree | `b2e6a423cb02acc64611d4e52dec047787345b85` / `94c524dfd5739feb648329f1154402f1b4844e30` |
| Delivered archive SHA-256 | `d9a7988509cd66b9c52b562debe151f7a4aa2f0a4676f559dd04fded1af213b8` |
| Diff SHA-256 | `3ef376898f50ca865586670592e76993524e28c7aa6697ea1161e9527c899958` |
| Background manifest SHA-256 | `a2b295083c73caf7d3d542ff24911967d0bf712609a5bd3f54ba53bf0683af20` |

The remote poster branch did not contain this commit. The supplied delivery bundle contains it with complete history and passes Git verification. All three attached reports match the archived copies byte for byte. Applying the archived diff to its base in a separate Git index yields the declared tree. All archive checksums pass. The final checkout is clean.

Twenty-six protected files were compared byte for byte with the base, including the frozen NREFT document, historical notebook, detector/rate/statistical source and Ge nuclear dataset/loader. No physical mass database, NREFT transformation or nuclear-data repair was introduced. Full hashes are in `SOURCE_IDENTITY.json`.

## Independent physics and numerical references

`bg_reference.py` imports no production physics. Expected beta spectra use 60-decimal mpmath Coulomb log-gamma functions and full-support normalization. Bateman activities and run averages use separate 80-decimal closed forms, including equal/near-equal half-lives and the actual stiff Ge68/Ga68 pair. Expected ES/CEvNS cross sections use explicit analytic kernels and total cross-section integrals, independently reconstructed unit conversion, a trigonometric Helm form factor and an independently assembled matter potential. Solar continuum references integrate the raw piecewise-linear flux density with the physical threshold clipped inside each source interval, rather than using production neutrino quadrature.

Tests cover scalar/vector behavior, full beta support, isotope/branch and shell probabilities, exact line-bin CDF probabilities, gamma/annihilation outcome conservation, cooldown versus acquisition means, explicit Ar activity/mass/transfer/veto factors, electron-shell occupancy, ER bypass of quenching, NR folding, observed residuals, component overrides, exposure scaling, immutable inputs and checksum failures. Agreement of a mathematical approximation with these references is not a measurement of its physical accuracy.

### Beta and EC source checks

[NNDC's evaluated Ar39 decay cards](https://www.nndc.bnl.gov/nudat3/getdecaydataset.jsp?dsid=39ar+bM+decay+%28268+y%29&nucleus=39K) confirm 7/2−→3/2+, first-forbidden unique, Q=565(5) keV and a ground-state branch of one. I downloaded the official BetaShape 2.4 package independently, verified SHA-256 `fa0a90094084b417c3852bcc93adcadba097d6d05187516c9aae1241595ad30e`, ran its official benchmark (62 file comparisons, no differences), and ran a separately transcribed evaluated input. No Coder regeneration/export implementation supplied the reference spectrum.

All **5,651** independently generated Ar39 energy/density pairs are identical to the portable production table. Its unnormalized trapezoidal integral is `0.9999946000463571`; the runtime correctly normalizes over the complete endpoint. The independently integrated piecewise-linear mean is **218.0114155 keV**, consistent with the evaluated 218.8±2.1 keV. Ar39 is not replaced by an allowed shape. Its low-energy calculation still inherits the tool's physical corrections, source endpoint uncertainty and original calculation-grid interpolation; requested 0.1-keV output is not an independent 0.1-keV physical calculation.

[DDEP H3](https://www.lnhb.fr/nuclides/H-3_tables.pdf) distinguishes atomic Q=18.591 keV from the 18.564-keV electron endpoint; the chosen endpoint and 12.312-year half-life are correct for the cited snapshot. The H3 approximation and positron signs/limits pass independent numerical checks. This does not certify omitted exchange, molecular, finite-size or radiative corrections.

The approximate H3 mean is `5.68693825` keV, also consistent with the cited DDEP mean5.68±0.01 keV.

The [LBNL binding-energy table](https://xdb.lbl.gov/Section1/Table_1-1a.htm) confirms the daughter K energies: Ga 10.367, Cu 8.979, Mn 6.539 and Zn 9.659 keV. Fully contained Fe55 EC therefore belongs at **6.539 keV**, not the Mn K-alpha photon near 5.9 keV. L1/M1 representative energies and Ge stepping energies also match that table.

Cached primary DDEP entries independently confirm dominant K/L/M probabilities, Ge/Ga and Fe half-lives, Zn EC48.35%/50.23%, Zn beta+1.421%, and dominant Ga beta+87.68%. M+ correctly includes reported M plus the higher-shell residual. Rounded Zn/Ga complete branch totals are explicitly retained and normalized; they are not unexplained normalization changes. Independent full-energy integration of each assembled decay reproduces its mean activity, including all gamma/annihilation outcomes.

**Source-validation limit:** direct DDEP PDF downloads failed with HTTP503/timeouts, although primary cached numerical entries were available. A complete branch-by-branch transcription audit, especially every weak Ga68 excited-state entry, was not possible. Original PDF checksums are not frozen in the snapshot, and Fe55's evaluation metadata says “current table”. Those gaps are not closed by source regeneration, branch-sum tests or agreement with the Coder's outputs. `PRIMARY_SOURCE_CHECKS.md` records checked versus unverified inputs explicitly. IAS direct source downloads returned403; cached primary Bahcall tables were spot-checked, while full frozen raw-table processing was independently verified.

### Time evolution, scattering and normalization

The chain uses activities correctly: `dAd/dt=b lambda_d Ap-lambda_d Ad`. Independent closed forms confirm the daughter growth term, equal-half-life limit and acquisition averages; no endpoint activity is substituted for a run mean. Initial Ga activity equal to initial Ge activity is a declared secular approximation followed by exact evolution, not an extra independently counted equilibrium component.

The absolute ES/CEvNS tests verify one conversion from GeV^-2 to cm², one GeV-to-keV Jacobian, target counts and seconds/day. With a synthetic unit form factor, the leading CEvNS integral is `GF² Qw² Eν²/(4pi)`; free-electron ES agrees with its analytic finite-endpoint integral. Survival probabilities reproduce independent vacuum, matter and high-energy limits. Raw pp/B8 source integrals and normalization agree; B16-GS98 totals agree with [Vinyoles et al. Table6](https://arxiv.org/html/1611.09867v4).

The AAr/UAr activity values independently match [DEAP](https://arxiv.org/abs/2302.14639) and [DarkSide](https://arxiv.org/pdf/1510.00702), including the mBq-to-Bq conversion. They are source-specific measurements, not LEGEND measurements. Full-source Ar intensity multiplies activity, 86400 s/day, Ar/Ge mass ratio, transfer and veto once each; efficiency then enters the assigned ER response once. Exposure scaling is distinct from rebuilding the run-averaged source activity.

## Independent poster comparison

A separate full-source-to-observed-window integrator (`compare_poster.py`) uses independent beta functions, closed activity averages, direct atomic-line Gaussian CDFs, independently coded neutrino kernels, raw solar spectra, independently coded Lindhard mapping and observed acceptance. It does not call the production response folding or plotting helper to generate expectations.

For 365.25 kg-day and 0.1–15 keV observed:

| Scenario | Independent total | Delivered example total |
| --- | ---: | ---: |
| Natural Ge / AAr | 2918.49896019 | 2918.49896180 |
| Enriched Ge / UAr | 110.81592835 | 110.81592811 |
| Poster, cooldown and assumed veto | **83.08333164** | **83.08333140** |

The poster independent components are H3 `62.38502223`, Ge68 `9.82103115`, Ga68 `0.87714657`, Zn65 `2.28995282`, Fe55 `5.52204519`, Ar39 `0.0184575375`, CEvNS `0.417681325`, ES `0.0384948281`, and observed flat `1.7135` events. Full 27-component comparisons are in CSV. Tritium and EC agree to floating-point precision; Ar39 differs by `7.49e-10` relative; solar terms differ by approximately `5.1e-7`/`7.3e-7` relative. These digits describe numerical agreement of the assumed model, not physical predictive precision.

The delivered refinement tool also runs successfully: its maximum component-integrated change is `7.313e-7` relative. Independent source quadrature refinement from order8 to12 changes ES by `3.38e-8` relative and CEvNS negligibly. Plots were regenerated and visually inspected; benchmark assumptions and observed-energy/year units are visible.

Doubling the independent recoil quadrature changes CEvNS by `1.23e-10` relative and ES by less than `4e-14`, supporting the stated numerical comparison.

## Reproduced implementation defects

**BG-AUDIT-001 — default Ar39 count integration fails.** A standard `ElectronRecoilResponse` with its default `IntegrationPolicy` cannot integrate the assembled Ar39 component over full support, even with `NoSmearing`; `expected_counts([-1,566], exposure)` raises `IntegrationWarning`, promoted to an exception by design. It also fails for the ordinary `[0.1,15]` window with the supplied Gaussian/Fano response. The tabulated shape has 5,651 knots, but `build_benchmark` supplies no table-knot truth breakpoints for Ar39. Direct QUADPACK diagnosis over full support gives relative actual error about `1.9e-8`, estimated error about `2.9e-6`, against requested `1e-8`. The source intensity and table normalization themselves are correct.

**BG-AUDIT-002 — default solar ES count integration fails.** Using the supplied Gaussian/Fano response and the poster configuration, `model.components['solar_ES'].expected_counts([0.1,15],365.25)` exhausts 300 quadrature subdivisions and raises. The fixed neutrino quadrature produces many recoil-endpoint kinks; the few supplied breakpoints do not make this default public fold reliable. Independent scalar kernels and the separate plotting fold pass. This is a count-assembly/numerical integration defect, not a wrong weak cross section.

Both failures are explicit exceptions, not silently accepted erroneous counts. They are scoped numerical integration defects; neither invalidates the independently reproduced plotting output. The repository tests use an explicitly looser `IntegrationPolicy(epsrel=1e-5,epsabs=1e-12)` for their factory count checks, which explains why 875 tests pass without exercising these default-policy failures. **The remedy must preserve accuracy requirements and supply a working public count path; tests were not relaxed and production was not repaired.**

## Regression results and acceptance boundary

Runtime: Python3.12.14, NumPy2.3.5, SciPy1.17.0, pytest9.1.1, mpmath1.4.1, Matplotlib3.10.8. Warnings were errors.

| Suite | Result |
| --- | --- |
| Unchanged repository | **875 passed**, 89.40 s |
| New independent benchmark suite | **179 passed, 4 failed**, 33.87 s; four manifestations of BG-AUDIT-001/002 |
| Unchanged earlier pre-rate/rate/analysis/background-SI suites | **652 passed, 1 failed**, 202.34 s |
| Unchanged independent Ge suite and diffraction-node checks | **573 passed, 3 failed**, 10.52 s |

The earlier architectural failure is the already reported duplicated frozen convention-ID literal. The three Ge failures are the already reported source J-closure deficit and Ge70/72 float64 diffraction-node negativity. These files are byte-identical to the base; they are inherited limitations, not introduced by this milestone. The unchanged detector's previously documented remote Gaussian-tail limitation was not re-audited.

**Acceptance:** use the supplied plotting route only as an explicitly illustrative benchmark, retaining activities/history, source-shaped LAr transfer, 1% veto, all-or-nothing containment, L1/M1 grouping, surrogate targets and uncalibrated low-energy response assumptions. Fix the two public-fold failures before claiming normal default `BenchmarkModel.expected_counts` integration is accepted. No quantitative LEGEND sensitivity or neutrino-floor acceptance is granted: measured composition/activity/transport/response are absent; Be7, pep, CNO and hep are omitted; source/model uncertainties are not propagated. Those missing physical inputs are distinct from the two reproducible implementation defects.
