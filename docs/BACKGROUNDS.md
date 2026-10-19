# Generic backgrounds and explicit response routing

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
not automatically applied again. This paragraph describes the original
observed-background models below; truth inputs use the explicit routing port
in the following section.

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

## NR truth, ER truth, and already observed inputs

`wimp_legend.energy.EnergyDomain` has three distinct enum members, `NR`, `ER`,
and `OBSERVED`, meaning keV_nr, keV_er, and keV_obs. Strings and booleans are
not interchangeable with these tags. The background module depends on a
structural response protocol, not a concrete detector implementation.

| Input | Public interface | Processing |
| --- | --- | --- |
| NR truth continuum/lines | `TruthBackground(EnergyDomain.NR, response, ...)` | Assigned response must declare `input_energy_domain=EnergyDomain.NR`; NR quenching, migration and efficiencies belong to that response |
| ER truth continuum/lines | `TruthBackground(EnergyDomain.ER, response, ...)` | Assigned response must declare `input_energy_domain=EnergyDomain.ER`; `ElectronRecoilResponse` maps E_er directly to its E_ee mean, with no NR quenching configuration |
| Already observed continuum | Existing `FlatBackground` / `CallableBackground` | Integrate observed rate directly; never refold or reapply efficiency |
| Already observed bin counts | Existing `BinnedBackground` | Use the declared exposure normalization directly |
| Already observed point lines | `ObservedLineBackground(lines, provenance=...)` | Assign observed energies directly to bins; no response or efficiency argument |

The common output remains `expected_counts(observed_bin_edges_keV,
exposure_kg_day)`: expected events in observed bins. `CompositeBackground`
accepts any mixture of these components. Each truth component owns its assigned
response; the response used for the signal is never applied to B again.

`TruthBackground` accepts `spectrum`, `spectrum_energy_range_keV`,
`spectrum_breakpoints_keV`, `lines`, and `provenance`. The spectrum is a scalar
callable returning events/(kg day keV_nr) or events/(kg day keV_er), according
to its explicit domain. The range/breakpoints apply only to this continuum;
every line has its own independently specified energy. A spectrum requires an
explicit finite nonnegative range with two increasing endpoints. Line-only
inputs do not take a continuum range. Negative/nonfinite rates are rejected.

`MonoenergeticLine(energy_keV, rate_per_kg_day, energy_domain)` records the
**integrated** line intensity in events/(kg day), not a differential density.
Truth energies are nonnegative; finite negative observed energies are allowed.
All lines must match their component's declared domain. For a truth line l,

\[
B_{jl}=X\,r_l\,\epsilon_{\rm truth}(E_l)
\int_{a_j}^{b_j}K(E_{\rm obs}\mid E_l)\epsilon_{\rm obs}(E_{\rm obs})dE_{\rm obs}.
\]

This calls `response.expected_line_counts(...)`, which uses bin probabilities.
It never approximates a delta function with ordinary truth quadrature. The
observed-efficiency integral inside a Gaussian bin may still use quadrature.
`DetectorResponse`, `ElectronRecoilResponse` and `KernelResponse` support this
line port. `MatrixResponse` supports NR continua, but its within-truth-bin
averages do not specify exact point-line migration: line routing is rejected.
An external ER kernel/matrix may implement the explicitly ER-tagged structural
port; it must supply its own line method if lines are requested.

Observed point lines use [a,b), with the final right edge included. A broadened
observed peak must be supplied as an observed continuum or histogram, not as an
unbroadened point line awaiting a second smearing step. Window/tail losses are
never normalized away. Input sequences are copied to tuples; callbacks and
external response objects must obey the existing pure-input contract.

```python
from wimp_legend.energy import EnergyDomain
from wimp_legend.backgrounds import TruthBackground, MonoenergeticLine, CompositeBackground, FlatBackground
from wimp_legend.detector import DetectorResponse, ElectronRecoilResponse, ConstantQuenching, NoSmearing
nr = TruthBackground(EnergyDomain.NR,
    DetectorResponse(ConstantQuenching(.5), NoSmearing()),
    spectrum=lambda E_nr: 1., spectrum_energy_range_keV=(0, 4))
er = TruthBackground(EnergyDomain.ER, ElectronRecoilResponse(NoSmearing()),
    lines=(MonoenergeticLine(1.5, 2., EnergyDomain.ER),))
background = CompositeBackground((nr, er, FlatBackground(.1)))
```

Every numerical value in this example is synthetic. This is response routing,
not a radioactive decay, neutrino, detector calibration or uncertainty model.
