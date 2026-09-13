# Particle-kernel API

The normative theory specification is [NREFT_CONVENTIONS.md](NREFT_CONVENTIONS.md).
This implementation covers its Wilson convention and Anand et al., PRC **89**,
065501 (2014), Eq. (38), only. The numerical rational factors in the kernel are
transcribed from that equation; no physical nuclear data are included.

## Public interface

Import from `wimp_legend.interactions`:

- `CANONICAL_OPERATORS`: immutable tuple `(1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)`.
- `NREFT_CONVENTION_ID`: the single implementation definition of the ID frozen
  in the normative document. The interactions layer owns this constant; no
  new shared layer or dependency on targets is needed.
- `WilsonCoefficients.from_isospin({operator: (c0, c1)})`.
- `WilsonCoefficients.from_proton_neutron({operator: (cp, cn)})`.
- `coefficients.get(operator, tau)` and `coefficients.to_proton_neutron()`.
- `RESPONSE_CHANNELS`: immutable tuple of the eight ASCII names below.
- `particle_response(channel, coefficients, *, tau, tau_prime, Q, V, j_chi)`.

`WilsonCoefficients(values)` also accepts canonical isospin pairs;
`WilsonCoefficients()` represents all-zero coefficients. The named factories
make input semantics explicit. The frozen dataclass copies all pairs into
tuples in a read-only `values` mapping. Omitted canonical operators become
zero. O2, unknown operators, noninteger indices, complex/nonfinite values,
booleans, and malformed pairs are rejected. Coefficients are real floating-point
values in GeV^-2, without weak-scale rescaling.

Conversion explicitly uses `c0=(cp+cn)/2`, `c1=(cp-cn)/2`, with inverse
`cp=c0+c1`, `cn=c0-c1`. The inverse returns a new read-only mapping for all
canonical operators. A conversion outside finite floating-point range fails.
`coefficients.convention_id` is fixed and cannot be supplied or replaced by
the caller. There are no automatic convention adapters.

## Response evaluation

The eight channels are `M`, `Sigma_double_prime`, `Sigma_prime`,
`Phi_double_prime`, `Phi_double_prime_M`, `Phi_tilde_prime`, `Delta`, and
`Delta_Sigma_prime`. Evaluation returns an ordered isospin entry, in GeV^-4.
The two mixed channels retain their ordered products and signs.

`Q=q^2/m_N^2` and `V=v_T_perp^2` are dimensionless independent inputs;
velocities are in units of c. Both must be real, finite, and nonnegative.
The kernel contains no numerical nucleon mass, target mass, or WIMP mass.
WIMP spin `j_chi` is required and must be a finite nonnegative integer or
half-integer: `2*j_chi` must be exactly integral, without tolerance or rounding.
Its spin factor must also fit in finite floating-point range.

NumPy broadcasting applies to Q and V, including constant-response channels.
A zero-dimensional output is a Python float; otherwise a fresh NumPy array
has the broadcast shape. Caller-owned inputs are never mutated. Overflow or
invalid response arithmetic raises an error instead of returning nonfinite
values. The evaluator does not include the external Q factors of Eq. (40),
nuclear spin averages, nuclear W functions, isotope weighting, or rate factors.

```python
from wimp_legend.interactions import WilsonCoefficients, particle_response

# Synthetic coefficients for illustrating the API, not benchmark input data.
c = WilsonCoefficients.from_proton_neutron({1: (2e-7, 2e-7)})
r = particle_response("M", c, tau=0, tau_prime=0,
                      Q=0.01, V=1e-6, j_chi=0.5)
# r = (2e-7)**2 GeV^-4
```

## Validation and limits

The full suite passes with `python -m pytest -q -W error`: **346 passed**, no
skips, failures, or warnings (Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0,
pytest 9.1.1). In the execution runtime, pytest was installed in a temporary
directory exposed through `PYTHONPATH`; no project dependency files changed.

The 155 new tests cover analytic p/n identities, all single-operator limits,
exact O1 behavior, interference signs and order, O12/O15 and O4/O6 cancellation,
spins 0, 1/2, and 1, kinematic boundaries, broadcasting, invalid inputs, and
immutability. A separate test-side sparse bilinear table expands Eq. (38)
independently of the production formulas and helpers, comparing all eight
channels and all four ordered isospin pairs on randomized coefficients and
momentum/velocity grids. An import-boundary test guards layer independence.
Existing nonempty tests were preserved; the empty interaction test skeleton
was populated. These checks validate the specified analytic implementation,
not the physical normalization of any nuclear dataset or absolute rate.

## Before the rate-layer run

Require exact nuclear-dataset convention-ID agreement and authoritative W
normalization/provenance. Implement the external Q factors for the five
Eq. (40) channels exactly once, along with the stated spin-average and
cross-section prefactors. Preserve ordered nuclear and particle entries.

Supply and record `m_N_reference_GeV` explicitly in analysis configuration;
changing it requires numerical regression checks but does not alone change
the analytic convention ID. Assemble Q and V consistently with elastic
kinematics and make rate-unit conversions and isotope counting explicit.
Then validate the absolute O1/SI benchmark, including the coherent zero-q
limit, before extending production rate calculations to other operators.

The historical Mathematica input-interface discrepancy remains isolated as
documented in the normative specification. No external-package comparison,
Wilson matching, nuclear contraction, cross section, or rate is implemented
in this run. No further discrepancy against the corrected specification was
found.
