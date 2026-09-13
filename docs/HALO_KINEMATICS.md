# Halo and elastic kinematics layer

Implemented against `refactor` at `0e6195befa55f7340612b48514b84c20189e2589`.
Scope: halo speed distributions, speed moments, elastic two-body kinematics,
and their tests. Nuclear responses, interactions, detector response,
backgrounds and statistics are unchanged. The empty `rates/kinematics.py`
placeholder is unchanged; new code should import package-level `kinematics`.

## API and dependency direction

- `halo.StandardHaloModel(v0, vesc, v_lab)` retains the existing constructor.
  All three arguments are scalar speeds divided by c. The model is frozen and
  copies scalar array inputs, preventing mutation through an external array.
- `StandardHaloModel.from_km_s(v0_km_s=220, vesc_km_s=544, v_lab_km_s=266)`
  is the explicit dimensional boundary and benchmark constructor.
- Model properties: `normalization_3d`, `z`, `max_lab_speed`, and
  `integration_breakpoints`. Methods: `galactic_speed_pdf(v)` and
  `lab_speed_pdf(v)`.
- `halo.velocity_moment(speed_pdf, v_max, v_min_c, *, n, breakpoints=(),
  epsrel=1e-9)` computes reference adaptive quadrature. The power is mandatory.
- `halo.build_velocity_integral(speed_pdf, v_max, *, power=-1,
  n_points=10000, breakpoints=())` retains the fast approximate callable API.
- `halo.mean_inverse_speed(speed_pdf, v_max, *, n_points=10000,
  breakpoints=())` returns the same table as `power=-1`.
- `kinematics.reduced_mass_GeV(m_chi_GeV, m_target_GeV)`.
- `kinematics.momentum_transfer_GeV(E_nr_keV, m_target_GeV)`.
- `kinematics.minimum_speed_c(E_nr_keV, m_chi_GeV, m_target_GeV)`.
- `kinematics.maximum_recoil_energy_keV(m_chi_GeV, m_target_GeV, v_max_c)`.

`shm` depends only on exact unit constants and generic numerics. `integrals`
accepts a PDF callable and support information; it imports no SHM or interaction
model. `kinematics` depends only on exact energy conversion and NumPy. Future
rate assembly can consume all three without halo/interaction coupling.

```python
from wimp_legend.halo import StandardHaloModel, velocity_moment, mean_inverse_speed
from wimp_legend.kinematics import minimum_speed_c

halo = StandardHaloModel.from_km_s()
vmin = minimum_speed_c([1., 5., 10.], m_chi_GeV=10., m_target_GeV=76*.9315)
eta_reference = velocity_moment(
    halo.lab_speed_pdf, halo.max_lab_speed, vmin, n=-1,
    breakpoints=halo.integration_breakpoints,
)
third_moment = velocity_moment(
    halo.lab_speed_pdf, halo.max_lab_speed, vmin, n=3,
    breakpoints=halo.integration_breakpoints,
)
eta_fast = mean_inverse_speed(
    halo.lab_speed_pdf, halo.max_lab_speed,
    breakpoints=halo.integration_breakpoints,
)(vmin)
```

## Physical conventions

The Galactic velocity density is proportional to `exp(-u²/v0²)` inside
`u < vesc` and zero outside. This is a hard-truncated Maxwellian, not a
lowered Maxwellian with a constant subtracted from its interior density.
`v0` is the Maxwellian scale; each Cartesian component has untruncated
standard deviation `v0/sqrt(2)`.

The truncation precedes the Galilean boost. Angular integration gives support
`max(0, v_lab-vesc) <= v < v_lab+vesc`, including boosts above the escape speed.
At zero boost the Galactic PDF is used exactly. The speed density is in
`d(beta)` with `beta=v/c`, and integrates to one. It contains no local DM
mass density. Negative speeds and infinities return zero; NaN is rejected.
PDF scalars return 0-d arrays; array shapes are preserved.

The integral convention is `I_n(a) = integral_a^vmax beta^n f1(beta) d(beta)`.
There is no hidden flux weight: eta is precisely `I_-1`. If a downstream API
uses speeds in km/s, its moments are `I_n(km/s) = C_KM_S**n * I_n(beta)`;
in particular divide eta by `299792.458` to express inverse speed in s/km.
Scalar/array threshold shape is preserved. Negative thresholds clamp to zero,
positive infinity gives zero, and NaN raises ValueError.

The model uses a fixed boost magnitude: 266 km/s is a benchmark snapshot,
not an annual-average prescription. No annual modulation, gravitational
focusing, anisotropic angular observables, or DM density is assumed here.
Speed moments are sufficient only after the relevant angular dependence has
been eliminated; directional observables need additional halo information.

Kinematics assumes a stationary free target and nonrelativistic, elastic,
two-body scattering. Mass inputs are rest energies `m*c²` in GeV; momentum
output means `q*c` in GeV (equivalently momentum in GeV/c); recoil input/output
is true nuclear-recoil keV, not electron-equivalent keV. `c=1` in the kernels,
and `1 keV = 1e-6 GeV` exactly. Kernels broadcast NumPy inputs without imposing
a target composition or halo cutoff. Finite positive masses and finite
nonnegative energies/speeds are required; endpoint speeds >=c are rejected.
`minimum_speed_c` returns the formula without clipping to halo support, even
for inaccessible recoil energies. The caller must remain in the NR regime.

## Numerical choices and limits

- The SHM normalization uses the regularized incomplete gamma function to
  avoid cancellation of the usual two-term erf expression at small vesc/v0.
- Full angular support uses `exprel`, and partial support uses `expm1`, to
  avoid cancellation at tiny boosts and near the endpoint. There is no
  arbitrary small-boost approximation to the PDF. Numerically unresolved
  quadrature breakpoints within 64 machine epsilons of an endpoint are omitted.
- Reference moments rescale the integration variable to `[0,1]`, factor out
  `vmax**n`, and use a global relative error budget across supplied breakpoints.
  QUADPACK warnings propagate as errors rather than producing apparently valid
  divergent integrals. Supply known kinks or narrow-feature boundaries: generic
  quadrature cannot discover every feature of an unknown callable.
- For SHM, `f1(v)=O(v²)` at zero, so real powers `n>-3` are integrable.
  `n<=-3` requires a strictly positive threshold in this API, conservatively
  even for a distribution with a gap near zero. More general PDFs must satisfy
  their own integrability condition. No infrared cutoff is inserted.
- The retained table is a composite-trapezoid/linear-interpolation approximation
  for powers >=-1 and PDFs with `f1(v)=O(v²)` at zero. Lower powers use the
  reference API. This explicit restriction replaces previously unrestricted
  but unreliable tabulation. The first grid cell is integrated, not discarded;
  the escape endpoint uses its interior PDF limit. Interior discontinuities
  should use reference quadrature. Tail-relative errors can be large arbitrarily
  close to vmax; use the reference calculation there and for validation.
- No new model hierarchy or interaction-specific halo objects are introduced.
  NumPy/SciPy runtime dependencies and a pytest extra are declared in pyproject.

## Independent validation

Run `python -m pytest -q -W error` from the repository root.

Tests cover all requested invariants, plus:

- independent angular integration of the 3D truncated Maxwellian;
- zero, tiny, equal-to-escape, and above-escape boosts;
- independent analytic piecewise eta and unboosted incomplete-gamma moments;
- `I_2(lab) = I_2(Galactic) + v_lab²`;
- a non-SHM polynomial PDF, fractional moments, and `dI_n/da = -a^n f1(a)`;
- tabulation convergence, invalid inputs, divergent moments, immutability,
  scalar/empty/multidimensional arrays and mass/energy broadcasting;
- random reduced-mass invariants and 527 random endpoint inversions;
- equal-mass energy transfer, heavy-WIMP limits and q=2*mu*v at the endpoint.

With target mass `76*0.9315 = 70.794 GeV` and maximum speed `810 km/s`:

| WIMP mass (GeV) | Maximum nuclear recoil (keV) |
|---:|---:|
| 5 | 4.49806195 |
| 10 | 15.83422672 |
| 20 | 50.15344515 |
| 40 | 134.72325829 |
| 100 | 354.33157886 |

Benchmark reference eta(0) is `1032.122085366866` in inverse beta units.
These are regression results for the specified inputs, not physical constants.

## Findings in the prior implementation

The existing refactor SHM is correct at the stated benchmark to ordinary
floating-point precision, but its above-escape-boost branch has the wrong
lower support bound. At 700 km/s it produces a minimum density of about
`-0.179609` on a 10001-point support grid. Direct subtraction of nearly equal
exponentials produces up to 17.35% relative error in a tested interior grid
at a boost of 1e-12 km/s. Nonfinite parameters were also accepted.

The previous generic integrator inserts a lower cutoff and accepts arbitrary
powers without checking integrability. For divergent powers that returns a
cutoff-dependent finite number; for near-singular integrable powers, a uniform
grid can be inaccurate. Its default inverse-speed benchmark is otherwise a
reasonable finite-grid approximation.

Notebook references below use zero-based JSON cell indices:

- Cell 4 `make_eta` returns zero below its first speed sample, including eta(0).
  Executing the selected notebook functions confirms eta(0)=0 versus the
  reference value above. The current notebook correctly includes the minus
  sign in reverse cumulative integration. The separate historical
  `WIMP_O1_spectrum_150626.py` export lacks that sign and plots absolute rates,
  masking the negative eta; do not conflate that export with the current notebook.
- Cell 9 boosted PDF divides by the boost with no zero-boost branch, suffers
  subtraction cancellation, has no nonnegative-speed guard, and assumes
  boost < escape. `zeros_like(v)` also truncates the density for integer inputs.
- Cell 33 uses `vesc+v0` for recoil endpoints despite its boosted PDF using
  `vesc+v_E`. With the stated benchmark this reduces Emax by a factor
  `(764/810)² = 0.889645`, about 11%. Later cell 40 uses the correct endpoint.
- Cells 33 and 40 read mutable global masses/speeds inside computations;
  redefining global parameters can mismatch an already-built eta table.
- Naming mixes `E` in GeV and `E_keV`, and `m_n` versus `m_N`. Cell 35's
  `GeV_to_keV=1e-6` is misleading as a bare unit-conversion name (the rate
  Jacobian has different semantics). This work does not audit or change rate
  normalization, nuclear responses, or detector physics.

The notebook and historical export are deliberately retained unchanged.

## Decisions before targets/nuclear responses

1. Specify isotope composition and whether fractions are by atom number or mass.
2. Select isotope masses and provenance; distinguish atomic and nuclear masses.
   `A*0.9315 GeV` is only a regression approximation, not an isotope database.
3. Agree on response momentum units, normalizations, and the definition of
   nuclear versus nucleon reduced masses at the interfaces.
4. Set the interaction-to-moment contract, including powers of c and any
   reference-speed normalization; keep density/flux factors in rate assembly.
5. Decide whether the first scientific result uses a fixed boost, annual
   averaging, or modulation, and set a moment accuracy budget before optimizing
   or caching tables for parameter scans.
