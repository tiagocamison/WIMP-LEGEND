# Generic reference NREFT rate layer

Theory authority: [NREFT_CONVENTIONS.md](NREFT_CONVENTIONS.md), Anand et al.
Eq. (40), Eq. (50), and Eqs. (52–53). The analytic convention is unchanged.
This work begins only after the Phase-A checkpoint; see
[CODE_AUDIT_001_RESOLUTION.md](CODE_AUDIT_001_RESOLUTION.md).

## Dimensional derivation (before implementing constants)

A natural-unit cross section has dimension GeV^-2. Restoring length units gives
`1 GeV^-2 = (hbar*c [GeV cm])^2 cm^2`. Consequently a differential cross section
`d sigma/d E_GeV` in GeV^-3 converts to cm^2/keV by multiplying by
`(hbar*c)^2 * dE_GeV/dE_keV = (hbar*c)^2 * 1e-6`.
This is an energy derivative Jacobian, not its inverse.

For f1 normalized in d(beta), the number flux is
`(rho_chi_GeV_cm3 / m_chi_GeV) * beta * c_cm_s`. Hence

```
dR/dE_keV = targets_per_kg * (rho_chi_GeV_cm3 / m_chi_GeV)
           * c_cm_s * seconds_per_day * (hbar*c)^2 * 1e-6
           * integral[f1(beta) * beta * d_sigma/d_E_GeV] d(beta).
```

The dimensions are `(targets/kg)*(1/cm^3)*(cm/s)*(s/day)*(cm^2/keV)`:
events/(kg day keV). There is exactly one physical flux factor c because beta
is dimensionless. No hbar/time conversion or extra powers of c occur.

Constants are derived from exact SI definitions, not a rounded particle-physics
conversion lookup: h = 6.62607015e-34 J s, elementary charge = 1.602176634e-19 C,
c = 299792458 m/s. Thus `GeV_to_J = 1e9 * elementary_charge`,
`hbar*c [GeV cm] = h*c_cm_s/(2*pi*GeV_to_J)`.
Sources: [BIPM defining constants](https://www.bipm.org/en/measurement-units/si-defining-constants)
and [SI Brochure](https://www.bipm.org/en/publications/si-brochure), Chapter 4,
Table 8 for day = 86400 s and the electronvolt; metric prefixes are exact.
The floating-point evaluations of these exact definitions incur ordinary
rounding. All these conversions live in `constants.py`.

## Public API and ownership

All exports below are from `wimp_legend.rates`:

```python
NREFTRateConfig(m_N_reference_GeV)

transition_probability(dataset, coefficients, *, q2_GeV2, V, j_chi, config)

differential_cross_section_GeV_minus3(
    E_nr_keV, beta, *, m_chi_GeV, dataset, coefficients, j_chi, config)

integrate_speed_flux(
    speed_pdf, differential_cross_section, *, beta_min, beta_max,
    breakpoints=(), epsrel=1e-8, epsabs=0.0)

differential_rate_per_kg_day_keV(
    E_nr_keV, *, m_chi_GeV, rho_chi_GeV_cm3, targets_per_kg,
    dataset, coefficients, j_chi, config, speed_pdf, beta_max,
    breakpoints=(), epsrel=1e-8, epsabs=0.0)
```

The frozen config contains exactly one finite strictly positive scalar mass in
GeV, required with no default. It stores an owned float. The caller must retain
this value with dataset identity, coefficients, j_chi, halo input, density and
target count as analysis provenance. It does not own numerical integration
policy. Changing the numerical mass does not change the analytic convention ID,
but must trigger numerical regression tests (O11's mass scaling is tested here).

| Quantity | Units/meaning |
| --- | --- |
| `q2_GeV2` | Momentum squared, GeV² |
| `V` | v_T_perp² in units of c²; finite nonnegative |
| `beta`, `beta_min`, `beta_max`, `breakpoints` | Dimensionless speed v/c |
| `E_nr_keV` | True nuclear-recoil keV, not electron-equivalent energy |
| `m_chi_GeV`, config mass, supplied isotope mass | Rest energy in GeV |
| `j_chi`, supplied nuclear `spin` J | Spin in units of hbar |
| `coefficients` | Canonical real Wilson coefficients in GeV^-2 |
| `transition_probability` output | P_tot, GeV^-4 |
| `differential_cross_section_GeV_minus3` output | d sigma/d E_GeV, GeV^-3 |
| Generic speed-flux output and `epsabs` | GeV^-3; includes beta, but not physical c |
| `speed_pdf` | Density in d(beta), normalized by the caller |
| `rho_chi_GeV_cm3` | Local DM mass density, GeV/cm³ |
| `targets_per_kg` | Supplied number of this isotope's nuclei per kg detector |
| Final rate | events/(kg day keV) |

Input masses, density, count, and spin are scalars. Density/count may be zero;
masses must be positive. q²/V arrays broadcast in the probability function; recoil/speed arrays broadcast in the cross section.
Recoil arrays are integrated point by point in the rate, preserving arbitrary
shape and empty arrays. Zero-dimensional outputs are Python floats. All inputs
are validated without mutation; returned arrays are owned. Negative/nonfinite
physical inputs, booleans, complex inputs and strings fail. Speeds must be below
1; physical use still requires the nonrelativistic regime. beta_min can exceed
beta_max, in which case the integral is zero.

## Dependency direction and convention gate

`rates/recoil.py` consumes the existing package-level kinematics, interaction
kernel, nuclear dataset interface, constants, NumPy and SciPy. The generic
integrator accepts scalar callables. The independent foundations do not import
rates or each other. The empty `rates/kinematics.py` is retained, not duplicated.

Before evaluation, require exact equality of the dataset's convention ID with
`NREFT_CONVENTION_ID`. Canonical datasets must additionally declare labels
`('0','1')` and `response_units='dimensionless'`; inconsistent declarations fail.
No descriptive string is parsed to infer a convention. These metadata gates do
not independently certify that imported physical coefficients actually obey it.

This initial reference API deliberately requires all 8 channels × 4 ordered
isospin entries, even for zero Wilson coefficients or zero-length input arrays.
A physically vanishing entry must be supplied as an explicit zero. Missing
entries raise KeyError, not zero or an inferred transpose. Extra dataset entries
are not contracted. The rate layer performs no automatic convention adapters,
proton/neutron transformations, abundance weighting, or source normalization.

## Q, transverse velocity, and Eq. (40)

Q = q²/m_N_reference² is calculated as `(sqrt(q²)/m_N_reference)**2` to avoid
squaring the dimensional reference mass before division. The particle kernel
still receives only Q, V and its explicit Wilson/spin inputs. Nuclear W evaluation
uses the dataset's public q² interface, preserving its coordinate mapping,
mathematical domain and enforced source-supported range.

All ordered terms are summed and multiplied by `4*pi/(2*J+1)`, with J from the
supplied isotope. One additional Q multiplies exactly:
`Phi_double_prime`, `Phi_double_prime_M`, `Phi_tilde_prime`, `Delta`, and
`Delta_Sigma_prime`. No external Q is applied to M, Sigma-double-prime, or
Sigma-prime. Signed terms and signed totals remain signed; arbitrary synthetic
or inconsistent datasets need not yield physical positive cross sections.
This kernel never hides such input problems with abs() or clipping.

For cross sections, validated package-level kinematics supplies q and beta_min.
Accessibility is checked first. For accessible speeds, use
`V=(beta-beta_min)*(beta+beta_min)`. Exactly beta=beta_min>0 gives V=0.
Below threshold the cross section is zero; inaccessible points do not trigger
nuclear evaluation outside its domain, although metadata/completeness are still
checked. `(E,beta)=(0,0)` has an undefined differential cross section and raises
ValueError. At E=0 the rate integral is allowed when convergent: adaptive interior
quadrature never evaluates the singular beta=0 endpoint. No hidden cutoff is used.

The cross section is precisely `m_T*P_tot/(2*pi*beta²)` per GeV recoil energy.
The relativistic amplitude factor `(4*m_chi*m_T)^2` is absent.

## Numerical integration and limits

`integrate_speed_flux` uses SciPy adaptive quadrature in x=beta/beta_max, with
scaled PDF beta_max*f1(beta_max*x). It numerically evaluates the supplied cross
section at each sample. Known speed breakpoints are filtered to the accessible
interval and scaled to x; duplicates are harmless. No moment reduction is used,
and a nonpolynomial cross-section test exercises this generic contract.

`epsrel` must lie between 50 machine epsilons and 1; `epsabs` is a nonnegative
absolute budget in the natural-unit speed-flux integral. IntegrationWarning is
raised as an error. PDF samples must be finite nonnegative scalars, and cross
section samples finite real scalars; source errors propagate. Supplied PDFs
must be normalized in d(beta) and supported within beta_max. The implementation
does not add a second normalization integral or silently normalize a PDF.
Callers must provide important breakpoints: adaptive quadrature cannot guarantee
discovery of arbitrarily narrow features. Signed near-cancellation may require
an explicit absolute integration budget. No global cache/state is introduced.

Finite float arithmetic is checked for overflow/nonfinite output. q² and Q
underflow are explicitly rejected. Arbitrarily extreme products, true sub-float
rates, or cancellations across many signed R*W terms are not arbitrary-precision
calculations. This reference path prioritizes transparent assembly; repeated W
and R evaluation per sample is intentionally not a production scan optimization.

`targets_per_kg` is supplied rather than inferred from bare nuclear mass:
nuclear mass fractions are not automatically measured bulk detector-mass
fractions. For a mixture, the caller supplies each isotope's appropriate target
count and sums isotope rates; no abundance is inserted into W or P_tot.

## Synthetic analytic validation

All new nuclear data live exclusively in `tests/test_rates.py`. Every fixture
specifies isotope properties, provenance, momentum mapping, conventions and all
32 responses including zeros. No production physical Ge data are introduced.

For f1(beta)=3 beta²/L³ and a=beta_min in [0,L], direct integration gives

```
eta(a) = 3/(2L) * (1-a²/L²)
B(a) = integral_a^L f1(beta)/beta * (beta²-a²) d beta
     = 3L/4 * (1-a²/L²)^2.
```

For W_M^(00)=1, J=1/2, j_chi=1/2, and sole coefficient c:
P_tot=2pi*R, so the speed-averaged cross section is m_T times
`c²*eta` for O1; `c²*Q*eta/4` for O11; and `c²*B/4` for O8.
Multiplying by the independently derived SI conversion and target/density
factors gives closed-form expected physical rates. These test oracles do not
call production particle/moment helpers. E=0, accessible nonzero recoils, and
inaccessible recoils are checked together in shaped arrays.

Other checks isolate each of the eight channels at multiple Q/V values,
including all five outer-Q channels; compare unequal ordered interference
entries with explicit scalar contractions; test every metadata/missing-entry
gate; and check source-coordinate and validated-range propagation. A synthetic
coherent W gives P_tot=(cp*Z+cn*N)^2 at q=0 and, after integrating the cross
section over recoil, sigma=mu²*(cp*Z+cn*N)^2/pi. This is an analytic convention
check, not validation of a real nuclear dataset.

The generic callback integral also reproduces
`integral_0^L f1(beta)*beta*exp(beta/L) d beta = 3L*(6-2e)` and an explicit
piecewise PDF with a breakpoint. An SHM O1 check uses an independent analytic
erf expression for eta(0). Additional tests cover SI constants via Decimal,
shapes, ownership, explicit density/count scaling, invalid inputs, singular
quadrature, and foundation dependency boundaries.

Validation: **492 repository tests passed**, including **109 new rate cases**;
**154/154 unchanged independent audit tests passed** after Phase B. No skips,
warnings or failures; commands use `python -m pytest -q -W error` with the same
temporary test dependency path documented in the Phase-A resolution. Runtime:
Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0, pytest 9.1.1.

## Before subsequent milestones

Before a historical O1 comparison, align isotope and nuclear mass, W form factor
and normalization, the sigma/c1 convention, explicit nucleon reference mass,
halo parameters/support, rho_chi, nuclear-versus-visible recoil units, and target
counts/kg. Do not change the notebook to force agreement.

Physical Ge import additionally requires versioned primary nuclear coefficients,
source normalization/isospin verification, oscillator-length/momentum mapping,
missing-versus-zero channel information, any validated range, and correct isotope
mass/spin provenance. This synthetic milestone certifies none of those data.

An optimized elastic I_-1/I_1 backend must reproduce this direct numerical
reference for all operators/interference terms and relevant thresholds before
replacing it in scans. Inelastic processes, light mediators, other DM species,
external-package adapters, detector response, backgrounds and sensitivity remain
future work. No frozen equation or convention ID was changed.
