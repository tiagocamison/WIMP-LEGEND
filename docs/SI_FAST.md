# Optional factorized SHM O1 calculation

The generic NREFT engine remains unchanged and is the default. Set the one
master-configuration field to opt in:

```json
"si_rate_engine": "factorized_shm"
```

This field belongs to `signal` in `analyses/legend_poster/config.json`. It can
be edited/saved in `notebooks/poster_setup.ipynb`; rerun the analysis notebook
or CLI normally. `"generic"` restores the reference implementation. Older
configurations lacking the field retain the generic behavior.

Only O1 is accelerated. O8/O11 and any mixed-operator call use the original
generic engine. The low-level gate tests every non-O1 stored coefficient for
**exact equality to zero**, without amplitude thresholds, spin/nuclear selection
rules or assumptions about isospin. All-zero coefficients are eligible too.

## Mathematics and conventions

Under the existing frozen half-isospin convention, define

\[
C_M(q^2)=\sum_{\tau,\tau'=0,1}c_1^\tau c_1^{\tau'}W_M^{\tau\tau'}(q^2).
\]

AFH PRC 89, 065501 (2014), Eqs. (40), (50), with the project's existing rate
assembly, then gives

\[
\frac{dR}{dE_{\rm nr}[\mathrm{keV}]}=
N_T\frac{\rho_\chi}{m_\chi}c\,(86400\ \mathrm{s/day})
(\hbar c)^2(10^{-6}\ \mathrm{GeV/keV})
\frac{2m_T}{2J+1}C_M(q^2)\eta(v_{\min}).
\]

All four ordered M entries are retained separately. Signed contributions and
the original AFH closure deficits are retained; there is no Helm approximation,
symmetrization, A^2 repair or clipping. `N_T` is the supplied single-isotope
nuclei/kg count, with composition assembled by the existing workbench. Masses
are GeV, coefficients GeV^-2, q^2 GeV^2, velocities beta=v/c, and the output is
events/(kg day keV_nr). The same core unit constants are reused.

For SHM eta, see Freese, Lisanti & Savage, Rev.Mod.Phys.85 (2013) 1561,
[arXiv:1209.3339v3 Appendix B Eq. (B4)](https://arxiv.org/html/1209.3339v3#A2),
with their truncation parameter beta=0. That parameter is distinct from our
speed variable beta=v/c. The implementation integrates the accepted local SHM
PDF and uses erf/erfc branches, a stable zero-boost expression, and 16-node
Gaussian evaluation of narrow exponential-difference integrals. This avoids
small-boost and endpoint subtraction cancellation without replacing the halo
by an unboosted approximation or inventing a cutoff. It uses no interpolated
eta table or global mutable physics state.

The standard SI proton conversion is unchanged:
`sigma_p_cm2 = GEV_MINUS2_TO_CM2 * mu_chi_p_GeV**2/pi * cp**2`.
The physical adapter converts sigma directly to coefficients before evaluating
the factorized rate; it does not evaluate another reference spectrum for this
conversion. Limits still obey `sigma_p_upper = mu_upper * sigma_p_ref_cm2`.
Pure-neutron and zero-coefficient models can use the coefficient dispatcher;
they cannot generally be represented by a finite neutron/proton ratio with a
zero reference proton cross section. This pre-existing distinction is preserved.

## Public interfaces

* `wimp_legend.halo.shm_eta.mean_inverse_speed_shm(halo, v_min_c)`:
  scalar/vectorized eta in beta units; negative thresholds clamp to zero,
  infinity is outside support, NaN raises, and the endpoint is zero.
* `wimp_legend.si_fast.is_o1_only(coefficients)`: exact eligibility predicate.
* `wimp_legend.si_fast.differential_rate_shm_per_kg_day_keV(E_nr_keV,
  halo=..., use_fast_si=False, **generic_arguments)`: the generic single-isotope
  arguments, except the explicit SHM owns the PDF, support and breakpoints.
  Opt-in plus exact O1-only coefficients selects factorization; otherwise the
  original generic integrator is called, with coefficients intact.
* `SIO1Signal(..., use_fast_si=True, halo=halo)`: physical SI adapter using the
  same conversion and cross-section limit method. Its supplied PDF/support/
  breakpoints must match that explicit SHM exactly. Custom PDFs remain supported
  by the default generic adapter; a contradictory fast request raises.

Both paths enforce the same canonical dataset convention, units, spin/config
and complete-entry checks before evaluation. Inaccessible recoils do not sample
nuclear functions. Fast evaluation samples the physically relevant M entries;
the generic engine still evaluates its full channel contraction. A required
nucleon reference mass remains explicit configuration, although O1 does not
depend on it analytically. No analytic convention ID or production rate formula
has changed.

## Caching and reproducibility

Truth rates remain keyed by isotope, halo, masses, couplings, engine and truth
grid/integration settings. They do not depend on abundance, exposure, background,
detector threshold or response. Folded rates are additionally keyed by detector
response and actual observed bin edges. A threshold on an existing edge can
reuse them; an off-grid threshold inserts a new edge and requires new integrals.
Background quadrature settings no longer invalidate signal folds. Detector
noise/efficiency changes require new folds but reuse truth evaluations. Changing
the SI option does not invalidate non-SI rates or folds. All cache keys still
include actual source/data/environment hashes; no old generic snapshot is
silently relabeled as a fast calculation.

The original standard snapshot is unchanged and remains an immutable generic
comparison. `tools/benchmark_si_fast.py --output DIR` validates 35 pure-isotope
folds and 84 SI limits against it, then measures matched generic/fast raw scans
and tracked allocation peaks. A changed physical master configuration fails
the benchmark's reference-configuration gate; it must not be compared to that
old snapshot as though assumptions were identical.

See `analyses/legend_poster/SI_SPEEDUP_REPORT.md` for measured performance,
validation tolerances, evidence and the independent tester handoff.
