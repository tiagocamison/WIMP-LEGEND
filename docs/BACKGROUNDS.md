# Generic observed-energy backgrounds

No radioactive, cosmogenic, neutrino, or LEGEND spectrum is supplied. The
existing radioactive/neutrino skeletons remain untouched. All empirical rates,
shapes, normalizations, uncertainties and provenance are external inputs.

The common `BackgroundModel.expected_counts(observed_bin_edges_keV,
exposure_kg_day)` returns an owned one-dimensional ndarray of expected events,
matching the detector count port. Differential rates are **events/(kg day
keV_obs)**. For exposure X in kg day:

\[
B_j=X\int_{a_j}^{b_j}b(E_{\rm obs})\,dE_{\rm obs}.
\]

Finite negative reconstructed-energy edges are permitted. Edges must increase
strictly; exposure is finite nonnegative. Counts and rate samples must be finite
nonnegative; booleans, complex, malformed shapes and nonfinite inputs fail.
Background rates already describe observed energy: detector efficiencies are
not automatically applied again.

## Public API (wimp_legend.backgrounds)

| Model | Contract |
| --- | --- |
| `ZeroBackground()` | zero counts after validating bins/exposure |
| `FlatBackground(rate_per_kg_day_keV_obs)` | counts = exposure × rate × bin width; explicit nonnegative rate |
| `CallableBackground(callback,policy=IntegrationPolicy(),breakpoints_keV_obs=(),provenance='user-supplied; unverified')` | scalar rate callback integrated independently in each bin; `.differential_rate(E_obs_keV)` supports scalar/array/empty inputs |
| `BinnedBackground(observed_bin_edges_keV,values,normalization,exposure_kg_day=None,provenance='user-supplied; unverified')` | exact-edge immutable binned input; explicit mode below |
| `CompositeBackground(components)` | tuple-owned independent components summed bin by bin; empty tuple means zero |

For `normalization='per_unit_exposure'`, values are **events/(kg day)** per bin,
and a fixed exposure is forbidden. The requested exposure multiplies them.
For `normalization='fixed_expected_counts'`, values are already **events**; an
explicit reference exposure is required. A request at any other exposure raises
ValueError; fixed values are never silently rescaled. Nonzero fixed counts at
zero reference exposure are rejected. Binned models require exact edge equality;
there is no implicit rebinning or interpolation.

Stored edges, values, breakpoints, and component containers are immutable tuples.
Returned arrays are independent copies. Custom component objects/callbacks must
remain pure. The default provenance string explicitly marks an unverified input;
replace it with a source/version/normalization description for real analyses.
These descriptive records do not certify any empirical dataset. Keep original
source data and uncertainty treatment separately with the analysis configuration.

Callable models use `IntegrationPolicy` from the generic numerical utilities.
IntegrationWarning is an error. Supply known narrow features/discontinuities in
`breakpoints_keV_obs`; adaptive quadrature does not promise to discover them.
Absolute tolerance has units events/(kg day), before exposure multiplication.
No signal, detector, target, or interaction module is imported by this layer.

```python
from wimp_legend.backgrounds import *
background = CompositeBackground([
    FlatBackground(.2),
    BinnedBackground([0,1,3], [.1,.4], 'per_unit_exposure',
                     provenance='synthetic example'),
])
B = background.expected_counts([0,1,3], 10.)  # [3.,8.] expected events
```

The example is synthetic; an observed histogram is not automatically a known
background template. Statistical inference here treats the supplied B as known;
empirical uncertainty requires a future/external nuisance-aware backend.
