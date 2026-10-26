# Poster analysis workbench — implementation report

## Identity and scope

Implementation commit: **`a819570f0f0669ab7b286c3dc77fd36a05cf77c4`**.
Implementation tree: **`07f365fc70426c6b938269864e4d34172fb3dc49`**.
This is byte-identical to the locally tested tree of commit
`150a83eb46b06fee5e6b40d61aa400cd23e5b354`; publication through GitHub assigned
the public commit identity above. The later delivery commit adds this report
and validation evidence plus hash-preserving Git attributes; it does not change
the executable implementation.

Only `poster-workbench` was advanced. The accepted core physics and all its
source datasets are unchanged. `poster`, `main`, `WIMP_O1_spectrum.ipynb` and
the pilot asset tables remain untouched. The old workbench notebook is preserved
byte-for-byte under `notebooks/historical/poster_workbench_pilot.ipynb`.

## Changes

* `analyses/legend_poster/config.json` is the one human-editable working input.
  It contains sourced isotope properties, explicit number/mass semantics,
  hypothetical enrichment, detector/background scenarios, halo and single-operator
  inputs, exposures, statistics, numerical profiles and plot selections.
* Small modules validate/save inputs, build existing package models, combine
  pure-isotope rates with correct per-kg weights, integrate detector bin
  probabilities, calculate limits, cache numerical arrays, and load/plot outputs.
  No background, nuclear, DM-response or statistical formula is duplicated.
* `poster_setup.ipynb` edits/validates/saves the master dictionary.
  `poster_workbench.ipynb` reads it, runs/reuses calculations, and inspects saved
  spectra, SI/O8/O11 sensitivities, exposures, counts, ROIs and convergence.
* Seven individual DirectDetectionPlots curves are frozen with raw Git/SHA-256
  hashes, upstream commit, references, 90% CL and observed/projection labels.
  They are separate from the signal calculation and explicitly may be outdated.
* Focused tests, a pinned optional notebook environment, CLI, usage guide and
  independent verification handoff accompany the implementation. `pyproject.toml`
  adds only the optional `workbench` extra.

The complete public entry points and scientific conventions are in `README.md`;
the independent-tester checklist is `TESTER_HANDOFF.md`.

## Scientific behavior

Number fractions are converted explicitly through `TargetComposition`.
Detector kilograms use neutral atomic bulk mass; recoil kinematics uses sourced
atomic masses converted to bare nuclear masses with an explicitly approximate
electron-binding correction. This keeps electrons out of the nuclear kinematics
while retaining their contribution to detector mass. Pure-isotope rates combine
to the configured mixture with no abundance or enrichment hard coding in core.

CEvNS sums the accepted weak kernel separately for each configured isotope.
Ge ES uses the exact total atom count and common Z/shell kernel. The standard
background factory and assigned ER/NR/observed routing are retained; discrete
lines use core line-bin probabilities. Already-observed residual/historical
backgrounds receive no extra response or efficiency.

O1 uses the existing physical SI adapter and `sigma_upper = mu_upper * sigma_ref`.
O8 and O11 are separate single-operator calculations with explicit spin,
proton/neutron coefficients and the frozen half-isospin convention. Their plotted
limits are `sqrt(mu_upper) * |c_p_ref|` in GeV^-2, never SI cross sections.
The nucleon reference mass is an explicit master input.

Default statistics are the existing exact Poisson-CLs construction in a
threshold-anchored 95% ROI, conditional on the configured observed window.
Optional full-window binned Asimov statistics remain separately labeled.
Historical/revised backgrounds share the new signal, exposure, ROI and statistics;
the historical observed detector selection is frozen, not recalibrated.

Every observed bin is reintegrated. PCHIP is restricted to adaptively refined
core truth rates; no coarse observed table is interpolated. Gaussian continuum
pruning has an explicit per-event probability bound; lines are unpruned.
Cache keys include relevant configuration, actual implementation/helper/data
hashes and numerical-library versions. Exposure or plotting changes reuse physics.
Manifests preserve resolved inputs, Git/dirty state, sources, numerical settings,
counts/statistical construction, convergence and numerical/figure checksums.

## Verification

Two non-overlapping final test runs cover **926 distinct tests**:

| Run | Command | Result |
| --- | --- | --- |
| Complete pre-existing suite, collected before the workbench test file was added | `python -m pytest -q -W error tests` | **889 passed in 2151.68 s** |
| Final focused workbench suite | `python -m pytest tests/test_poster_workbench.py -q -W error` | **37 passed in 57.83 s** |

This is not a claim of a single 926-test invocation. The final checkout can run
all tests together with `python -m pytest -q -W error`. No existing test or
assertion was edited or deselected in the completed pre-existing run. Logs are
included under `validation/`. Additional overlapping/superseded attempts are not
added to these totals.

The focused run uses the pinned fresh Python 3.12 environment. It executes every
cell of both actual notebooks in separate fresh IPython subprocesses, retains
rich outputs, verifies every cell executed, and checks that an unedited setup
round trip preserves the temporary master file. Managed CI blocks Jupyter NIC
discovery and ZMQ TCP/IPC services; no claim is made that a Jupyter server/kernel
transport passed here. Notebook mathematics and plotting are not mocked.

Tests cover configuration ownership/round trips/failures; isotope number/mass/
atomic-kg identities; ER/NR/observed routing and individual response overrides;
exact line probabilities; rebinning of independently integrated counts; real
0.005-versus-0.01-keV physical signal bins; SI/coefficient normalization and
reference scaling; mass-dependent limits; exact counting/optional binned methods;
exposure/cache reuse with expensive calls forbidden; cache corruption; frozen
comparison blobs; historical file identity; and saved-output/figure consistency.

The standard profile is exercised separately over all configured masses,
operators, isotopes, scenarios and exposures. Its saved provisional numerical
outputs, manifest, convergence summary and figures are described in
`validation/README.md`. They are a new reproducible benchmark, not promoted pilot
numbers. Plot layout and threshold structure were visually inspected.

The standard run produced **252 successful limits**, **105 isotope/signal checks**
and **27 background checks**. Maximum signal total refinement change was
`9.834080468537064e-6`; maximum significant-bin change was `0.02786773142963036`
(2.79%, within the configured 3%). Background checks agreed to numerical roundoff
for the accepted prepared kernels; this is not a physical-uncertainty estimate.
The published-commit rerun recorded **458 cache hits and zero misses**. Seven
figures are included in both vector PDF and PNG form. The optional binned backend
is tested separately and is not included in these standard counting limits.

One newly authored mass-dependence assertion initially demanded an unjustified
1% difference between two arbitrary masses. It was corrected to distinguish
different values without excluding a real sensitivity plateau. Existing tests
were unchanged. The initial notebook transport checks failed due to CI socket
restrictions; the fresh-process execution above replaced that launcher, without
removing cell-execution or numerical assertions.

## Remaining assumptions

The enrichment, activity history, detector calibration, LAr transfer/veto,
containment and exposures remain explicit provisional scenarios. Solar coverage
is pp/B8 only. Mass/binding and physical source uncertainties are not propagated.
AFH source closure deficits and unknown physical validated range are unchanged;
negative core rates raise rather than being clipped. Fixed halo and nuisance-free
statistics are benchmark choices. A denser mass grid and independent convergence
review are needed for final poster curves. External experimental assumptions are
not translated into our halo/statistical conventions.

Earlier core-engine audit limitations remain documented in the existing reports;
this analysis-only change neither repairs nor certifies them. In particular,
successful internal tests are not an independent nuclear-data or LEGEND detector
validation. The handoff lists checks for a separate scientific reviewer.
