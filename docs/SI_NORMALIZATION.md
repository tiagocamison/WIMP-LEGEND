# Physical SI normalization without a special rate engine

`wimp_legend.si` is a contact O1 adapter around the existing generic NREFT rate.
It does not alter `rates`, `statistics`, `sensitivity`, the operator kernels,
`GEV_MINUS2_TO_CM2`, or the frozen `NREFT_CONVENTIONS.md`.

The source is Anand, Fitzpatrick & Haxton, PRC 89, 065501 (2014),
[doi:10.1103/PhysRevC.89.065501](https://doi.org/10.1103/PhysRevC.89.065501),
Eqs. (15–18) and (50), as fixed in Section 17 of the normative document.
For a free proton, P_tot=(c1_p)^2 and integration over the elastic recoil
interval gives the normalization below. The existing conversion constant is
derived from exact SI definitions; no second numerical conversion is added.

## Scalar conversions

\[
\sigma_p[\mathrm{cm}^2]=C\,\frac{\mu_{\chi p}^2}{\pi}(c_1^p)^2,
\qquad
c_1^p=s_p\frac{\sqrt{\pi}\sqrt{\sigma_p}}{\sqrt C\,\mu_{\chi p}},
\quad C=\texttt{GEV\_MINUS2\_TO\_CM2},\quad s_p\in\{-1,+1\}.
\]

| Public API | Meaning |
| --- | --- |
| `si_proton_coefficient_GeV_minus2(sigma_p_cm2, *, m_chi_GeV, m_proton_GeV, proton_sign)` | Select the explicitly specified global sign of the dimensionful proton O1 coefficient |
| `si_proton_cross_section_cm2(c1_p_GeV_minus2, *, m_chi_GeV, m_proton_GeV)` | Invert a signed proton coefficient into its nonnegative physical cross section |
| `si_wilson_coefficients(sigma_p_cm2, *, m_chi_GeV, m_proton_GeV, fn_over_fp, proton_sign)` | Set only O1, with cn=(fn_over_fp) cp, through `WilsonCoefficients.from_proton_neutron` |

`fn_over_fp` **means c1_n/c1_p**. It is required: 1 gives equal p/n coupling,
0 gives proton-only coupling, and a negative ratio retains its relative sign.
`proton_sign` is also required. In an O1-only signal a simultaneous p/n sign
flip leaves rates unchanged; no relative sign is inferred from a cross section.
The canonical coefficients are c0=(cp+cn)/2, c1=(cp-cn)/2 in GeV^-2.
There is no factor 4, weak-scale coefficient or alternative full-sum convention.

Masses are required positive real scalars in GeV; no physical proton mass is
hard-coded. That mass is distinct from the required `NREFTRateConfig` nucleon
reference mass defining Q. The scalar conversions accept sigma=0 or cp=0;
a finite cn/cp ratio then sets both coefficients to zero. A pure-neutron model
with cp=0 and cn!=0 requires the generic Wilson API rather than this proton
reference parametrization. Nonfinite, malformed or unrepresentable nonzero
results fail explicitly. Representable subnormal outputs retain float precision.

## Reference signal and limit

`SIO1Signal` is a frozen callable configured with keyword arguments:

```text
sigma_p_ref_cm2, m_chi_GeV, m_proton_GeV, fn_over_fp, proton_sign,
dataset, rho_chi_GeV_cm3, targets_per_kg, j_chi, config,
speed_pdf, beta_max, breakpoints=(), epsrel=1e-8, epsabs=0
```

The reference cross section must be strictly positive. `signal(E_nr_keV)`
delegates to `differential_rate_per_kg_day_keV` and supplies only the constructed
O1 Wilson coefficients. Scalar/array/empty inputs retain the generic shapes.
All density, target-count, spin, halo and quadrature inputs remain explicit.
No isotope fraction, mass-to-count conversion or nuclear coefficient is added.
The dataset must satisfy all existing canonical ordered-entry/convention gates,
including explicit known-zero entries for channels unused by O1. Constructor
validation uses the generic empty-array path without sampling the PDF or W.

Pass the callable to the existing `SensitivityAnalysis`:

```python
result = analysis.observed_limit(signal, bins, exposure_kg_day,
                                 truth_energy_range_keV, observed_counts)
sigma_p_upper_cm2 = signal.cross_section_limit_cm2(result.limit.mu_upper)
```

At the default confidence level .90,

\[
\sigma_{p,90}=\mu_{90}\sigma_{p,\mathrm{ref}}.
\]

For fixed masses, ratio, dataset and other inputs, the template scales linearly
with sigma_ref, while coefficients scale with its square root. Changing the
arbitrary reference normalization by a positive factor leaves the final physical
bound unchanged within the chosen root/integration tolerance. This is a pure
cross-section rescaling; it does not reinterpret the statistical construction.
Retain the original LimitResult, especially its confidence level, method and
asymptotic/coverage notes. `mu_upper=+inf` preserves no_sensitivity as +inf.

The adapter describes one isotope. Sum separately configured isotope callables
with the same particle model/reference to build a mixture; the generic target
handling is unchanged. All regression datasets are synthetic coherent charges
or zero entries, with synthetic masses. No Ge nuclear response, decay chain,
experimental calibration or empirical isotope mass is included.
