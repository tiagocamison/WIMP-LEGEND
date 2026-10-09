# Generic detector response

No calibration, LEGEND enrichment, or physical nuclear dataset is supplied.
Lindhard quenching and a Gaussian resolution are configurable phenomenological
fallbacks. They are not an official LEGEND response. A collaboration may replace
any protocol with its own pure implementation without subclassing these models.

## Energy and count contract

The input callable takes **E_nr_keV**, true nuclear recoil energy, and returns
**events/(kg day keV_nr)**. Quenching maps it to a mean **E_ee_keV**, and
resolution maps that mean to **E_obs_keV**, reconstructed energy. Bin counts are
expected events, not a differential rate. Exposure is **kg day**.

For exposure X, the normative fold is

\[
S_j=X\int_{E_{\rm nr,lo}}^{E_{\rm nr,hi}}r(E_{\rm nr})
\epsilon_{\rm truth}(E_{\rm nr})
\int_{a_j}^{b_j}K(E_{\rm obs}\mid E_{\rm nr})
\epsilon_{\rm obs}(E_{\rm obs})\,dE_{\rm obs}\,dE_{\rm nr}.
\]

Truth energies must be finite nonnegative. Reconstructed bins may extend below
zero: a Gaussian on the real line produces negative reconstructed outcomes.
Bounds and bin edges are finite and strictly increasing. Truth range must be
supplied explicitly; it is never inferred from observed bins. Events outside
that truth range are excluded by the caller's model, and events outside the
observed window remain lost. Neither loss is renormalized away.

## Public API (wimp_legend.detector)

| API | Inputs and meaning |
| --- | --- |
| `QuenchingModel.mean_energy_keV(E_nr_keV)` | scalar/array true keV -> same-shape mean electron-equivalent keV |
| `ConstantQuenching(q)` | explicit 0<q<=1, E_ee=q E_nr; zero quenching unsupported |
| `LindhardQuenching(Z,k)` | explicit positive integer Z and positive k; `quenching_factor(E_nr_keV)` also available |
| `ResolutionModel.bin_probabilities(E_ee_mean_keV, observed_bin_edges_keV, *, observed_efficiency=None, policy=IntegrationPolicy(), observed_breakpoints=())` | shape mean.shape+(number of bins,), probabilities with optional observed efficiency integrated |
| `NoSmearing()` | deterministic reconstruction at E_ee |
| `GaussianResolution(sigma_keV)` | finite nonnegative constant keV or scalar callable of **E_ee_mean_keV**; `sigma_at(...)` preserves shape |
| `Efficiency.__call__(energy_keV)` | probability in [0,1], same scalar/array shape |
| `ConstantEfficiency(epsilon)` | explicit [0,1] |
| `CallableEfficiency(callback)` | scalar callback evaluated elementwise, supporting scalar or array callers |
| `DetectorResponse(quenching,resolution,truth_efficiency=None,observed_efficiency=None,policy=...,observed_breakpoints=())` | composed conditional response, no hidden quenching/resolution values |
| `DetectorResponse.bin_probabilities(E_nr_keV, observed_bin_edges_keV)` | includes observed efficiency; truth efficiency is applied only during truth folding |
| `CountResponse.expected_counts(truth_spectrum, observed_bin_edges_keV, exposure_kg_day, truth_energy_range_keV, *, truth_breakpoints=())` | common response port; returns owned one-dimensional expected-count ndarray |
| `KernelResponse(kernel,includes_efficiency,truth_efficiency=None,policy=...)` | bin-integrated custom kernel port, callback(E_nr_keV, observed_edges) -> probability vector |
| `MatrixResponse(truth_bin_edges_keV,observed_bin_edges_keV,matrix,includes_efficiency,within_bin_averaging,truth_efficiency=None,policy=...)` | fixed-edge matrix port; same expected_counts method |
| `MatrixResponse.apply_truth_counts(truth_counts)` | already integrated truth-bin counts -> observed counts |
| `IntegrationPolicy(epsrel=1e-8,epsabs=0,limit=300).integrate(callback,lower,upper,breakpoints=())` | finite positive-integrand adaptive quadrature; generic numerical utility |

None means unit efficiency at that stage, not a calibration assumption. No
input probability is silently clipped. Efficiency values have no tolerance:
0<=epsilon<=1 is enforced. Callbacks are scalar-sampled; returned booleans,
complex, nonscalar, negative or nonfinite values are rejected. Energy-array
inputs including empty and zero-dimensional arrays are supported. Configuration
numbers are owned floats/ints; array-like stored values become immutable tuples.
Callable closures and external models remain caller-owned pure contracts: their
captured mutable state cannot be deep-frozen by this interface.

## Mathematical models

The Lindhard fallback uses the formula supplied for GENERIC-ANALYSIS-001:

\[
\epsilon=11.5 E_{\rm nr,keV} Z^{-7/3},\qquad
 g=3\epsilon^{0.15}+0.7\epsilon^{0.6}+\epsilon,\qquad
 Q=\frac{kg}{1+kg},\qquad E_{\rm ee}=QE_{\rm nr}.
\]

The decimal coefficients define this phenomenological parametrization; no
isotope-dependent k or empirical numbers are supplied. Zero recoil gives zero
Q and mean. Positive Q is evaluated as `1/(1+(1/g)/k)` to avoid overflow of kg.

Gaussian migration is a normal distribution with mean E_ee and sigma evaluated
at that mean. Each bin probability is a normal-CDF difference. Bins in the
positive standardized tail use survival-CDF differences to avoid cancellation
near one. sigma=0 uses deterministic bins **[a,b)**; the final bin includes its
right edge. This convention matters for point-like truth inputs; it has measure
zero for regular differential spectra. No half-event splitting is used at a
zero-width boundary.

An arbitrary observed efficiency is integrated inside the bin after changing
to its conditional normal-probability coordinate. It is not sampled only at a
bin center. Constant efficiencies use the exact probability scaling shortcut.
Known observed-efficiency jumps should be supplied in `observed_breakpoints`.
True efficiency is sampled at E_nr, observed efficiency at E_obs.

## Custom ports and ownership

A kernel callback supplies the complete probability of reconstruction in each
observed bin. `includes_efficiency` is a **required bool**. If true, adding a
truth efficiency is rejected. If false, an optional truth efficiency may be
applied separately. A separate observed efficiency is not accepted by a
bin-integrated custom port: the callback must integrate it itself and declare
ownership. `includes_efficiency=True` conservatively owns both stages.

Matrix orientation is **M[j,i]**: observed bin j, truth bin i,
`observed_counts = M @ truth_counts`. Entries are nonnegative and <=1;
column sums must be <=1 within **8 machine epsilons** of probability roundoff.
The same sum tolerance applies to custom kernel vectors; no renormalization is
performed. Matrices are stored as immutable nested tuples. Edges must match
exactly; no implicit rebinning. `expected_counts` requires the full matrix truth
range to avoid an unstated partial-bin convention. The required nonempty
`within_bin_averaging` text records how each column was obtained. A matrix with
uniform-within-bin averages is not exact for arbitrary varying sub-bin spectra.
Extra observed efficiency is unavailable for the same sub-bin-information reason.
`apply_truth_counts` receives already integrated counts, including any separately
owned truth acceptance; the caller must honor this count-level contract.

## Numerics and example

Both truth folding and observed-efficiency folding use adaptive scalar SciPy
quadrature. `epsabs` has the units of each integral (before exposure for truth
folding; dimensionless conditional efficiency average for observed folding).
The returned count has events units. Relative tolerance must lie between 50
machine epsilons and 1; limit is a positive subdivision budget.
`IntegrationWarning` becomes an error. Nonfinite arithmetic fails explicitly.
Known breakpoints are copied, sorted/deduplicated, and filtered to each interval.
No guarantee is made that adaptive integration discovers arbitrary narrow
features in user callbacks. Supply their breakpoints and verify convergence.

For the two built-in monotone quenchers, deterministic bin crossings are split
explicitly. Constant Gaussian widths also insert truth crossings at each
observed edge and its +/-1,2,4,8 sigma positions. These are quadrature hints,
**not Gaussian truncation**. Callable widths or custom kernels must supply
additional truth breakpoints themselves. Repeated bin-by-bin folding is a
transparent reference implementation, not an optimized response cache.
Very extreme finite scales can still underflow intermediate integrands or fail;
this layer is not arbitrary precision and does not extend the rate product fix.

```python
from wimp_legend.detector import *
response = DetectorResponse(ConstantQuenching(.2), GaussianResolution(.1),
                            observed_efficiency=ConstantEfficiency(.4))
S = response.expected_counts(lambda E_nr_keV: 2., [0,.2,.6,1.2],
                             exposure_kg_day=100., truth_energy_range_keV=(0,5))
```

All numbers in this example are synthetic choices, not experimental data.
See `examples/generic_analysis.py` for the complete count-to-limit chain.
