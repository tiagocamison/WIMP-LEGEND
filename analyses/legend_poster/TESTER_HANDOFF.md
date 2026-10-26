# Independent verification handoff — poster workbench

Use a fresh editable checkout on `poster-workbench` and the environment described
in `README.md`. This document specifies checks, not an independent tester's
approval. The implementation author must not represent these checks as a second
person's audit.

## Entry points

* `analyses.legend_poster`: `load_config`, `save_config`, `validate_config`,
  `resolve_profile`, `run_analysis`, `load_run`.
* `models.target`, `isotope_weights`, `responses`, `build_background`.
* `analysis.observed_edges`, `fold_component`, `signal_model`, `pure_signal_counts`,
  `roi_slice`, `count_limit`.
* `plotting.comparison_curves`, `limit_series`, `plot_run`.
* `notebooks/poster_setup.ipynb` and `notebooks/poster_workbench.ipynb`.

The setup notebook owns editing/saving; the analysis notebook reads the saved
file. Both allow temporary configuration/output paths through environment variables.
Run `python -m pytest tests/test_poster_workbench.py -q -W error` for the focused
checks, including clean-IPython-subprocess execution of both notebooks. Set
`WORKBENCH_TEST_CACHE=analyses/legend_poster/cache` to reuse checksummed work.

## Critical scientific properties to challenge

1. For every isotope, verify
   `bulk_mass_fraction * pure_targets_per_kg == targets_per_kg`. Re-express a
   composition as nuclear mass fractions; final number weights/rates must agree.
   Natural includes five isotopes; the provisional enrichment has only 74/76.
   Check the NIST atomic masses, binding correction and kg normalization
   independently. The fraction meaning is deliberately not silently changed.
2. Compare a pure-isotope recoil rate at several energies directly with core
   `differential_rate_per_kg_day_keV` / `SIO1Signal`. Then sum pure saved bin
   integrals with atomic bulk weights and compare to saved mixtures. AFH data
   hashes, half-isospin convention and original closure deficits must be unchanged.
3. Check O1 sigma-to-coefficient round trips with explicit proton mass and
   neutron ratio. Rescaling sigma_ref by k must rescale reference signal by k
   and mu_upper by 1/k, leaving the physical sigma limit invariant. For O8/O11,
   rescaling coefficient_ref by k scales signal by k^2; coefficient limits must
   be invariant after the sqrt(mu) conversion. Repeat for a non-unity neutron
   ratio and alternative spin. Do not interpret coefficient limits as SI limits.
4. For H3/EC/argon/solar components verify assigned ER/NR response and exactly
   one efficiency. Observed residuals/historical curves must bypass response.
   Compare a line to the core CDF bin probability, not a sampled delta function.
   Change only one component's response preset and check isolation.
5. CEvNS must equal a sum of isotope-specific *existing* weak kernels with
   explicit nuclei/kg. ES must use the sum of nuclei/kg and common Ge shells.
   Compare the explicitly selected effective-isotope option separately.
6. Halve near-threshold bins (to 0.005 keV), refine truth grids and quadrature,
   and compare integrated rates, ROI boundaries and limits. A 95% bin-edge ROI
   is discrete; the exact expected-count Poisson construction can also have
   integer pseudo-observation transitions. Do not mistake these for a continuous
   asymptotic method. Failed configured convergence must raise.
7. Double exposure with the same activity history: signal/background counts
   double, spectra do not; physical limits are not generally proportional to
   1/exposure in nonzero background. Compare the exact backend call to saved
   mu/counts. Optional binned Asimov curves must remain distinctly labeled.
8. Change only plotting styles/exposures and verify no new expensive rate calls.
   Mutate a cache file and require checksum failure. Verify manifest file hashes,
   actual implementation hashes and Git state; run from the intended editable
   package, not an unrelated globally installed checkout.
9. Recreate figures from saved arrays without a physics run. Plotted ordinates
   must agree with `limits.csv` and component counts divided by actual bin width.
   Check vector PDFs, near-threshold line structure, comparison source labels,
   CLs and observed/projection distinctions. No pilot table may feed revised plots.
10. Execute both notebooks from fresh IPython subprocesses, including setup save/round trip
    using a temporary configuration; compare the saved config bytes before/after
    an unedited run. Preserve original pilot notebook/assets and WIMP notebook.

## Remaining physical assumptions

The 92/8 enrichment, detector calibration/resolution/efficiency, activities,
containment, LAr transfer/veto and exposure scenarios are provisional. Neither
the activity normalization nor a LEGEND geometry/transport model is validated.
Solar pp/B8 coverage is incomplete. Nuclear source functional forms have no
invented physical validity cutoff. Mass/binding and source uncertainties are not
propagated. Fixed halo and nuisance-free statistics remain benchmark choices.
The old observed background cannot be recalibrated by the new detector model.
External limit snapshots are not guaranteed latest-2026 results and use their
own assumptions. Smoke tolerances are for installation checks only. Final poster
plots need a denser mass grid and an independent numerical convergence study.

Previously documented core-engine limitations/audit findings are not repaired
or silently waived by this analysis-only change. Inspect the existing audit
reports before treating a successful workbench run as validation of all core
physics. This work adds no calibration or nuclear-structure evidence.
