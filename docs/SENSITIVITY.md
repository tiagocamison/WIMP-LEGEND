# Generic sensitivity orchestration

The driver delegates scientific models and statistical methods. It contains no
NREFT response, isotope handling, matching, detector calibration or background
physics. A truth provider may be any pure callable returning
**events/(kg day keV_nr)**, including an explicit adapter to the existing rate
function. No physical Ge data or LEGEND predictions are included.

## Public API (wimp_legend.sensitivity)

| API | Contract |
| --- | --- |
| `SensitivityAnalysis(response,background,backend)` | immutable holder of CountResponse, BackgroundModel, StatisticsBackend ports |
| `.observed_limit(truth_spectrum,observed_bin_edges_keV,exposure_kg_day,truth_energy_range_keV,observed_counts,*,confidence_level=.90,truth_breakpoints=())` | fold supplied unit-strength signal and background, validate actual integer n, delegate observed limit |
| `.expected_limit(truth_spectrum,observed_bin_edges_keV,exposure_kg_day,truth_energy_range_keV,*,confidence_level=.90,truth_breakpoints=())` | delegate backend-defined background-only expectation |
| `.expected_grid(truth_spectra,observed_bin_edges_keV,exposures_kg_day,truth_energy_range_keV,*,confidence_level=.90,truth_breakpoints=())` | named mapping of callable templates crossed with exposure vector -> tuple of GridResult |
| `AnalysisResult` | immutable signal/background count tuples, edges, exposure, truth range, and LimitResult |
| `GridResult` | immutable template_id, exposure_kg_day, analysis |
| `scan_parameter(signal_counts,*,parameter_domain,parameter_grid,background_counts,observed_counts,backend,confidence_level=.90)` | explicit nonlinear absolute-count scan; returns ParameterScan |
| `ParameterScan` | owned tuples of domain, parameters, margins, accepted flags; method, construction, approximation, confidence, status and notes |
| `.accepted_sample_ranges` | connected ranges of **accepted supplied samples**, not interpolated confidence intervals |

The default linear path is **S_i(mu)=mu*S_i(1)**. It is a strength limit on a
supplied template. It is not automatically a Wilson-coefficient/cross-section
limit: physical parameters may interfere or affect spectral shapes nonlinearly.
Every exposure is explicit in kg day and every integration truth range is
explicit in keV_nr. Fixed-count backgrounds still enforce their reference
exposure. A zero signal/exposure produces a no-sensitivity result, not a finite
sensitivity claim. Custom count providers and statistics need no subclassing.

## Nonlinear signal models

`signal_counts(theta)` returns absolute expected observed-bin events. The
explicit domain may include negative parameter values. The finite, strictly
increasing grid must span that domain exactly; it defines the resolution of the
scan. For each theta the backend tests the supplied signal at strength one.
Its `acceptance_margin` must return a finite scalar, nonnegative if accepted,
and `acceptance_metadata()` returns the method ID, test definition and
approximation. These are preserved in ParameterScan. Actual observations remain
integer counts and vectors must match.

The scan makes no monotonicity assumption. It reports all accepted samples,
including disconnected regions, and labels its status `completed_explicit_grid`.
It does not call an arbitrary crossing a unique upper limit, and cannot rule
out narrower allowed/excluded regions between samples. Refinement or a dedicated
parameter-specific statistical construction is needed for continuous confidence
regions. The built-in shape test profiles strength along each point's signal
shape; it does not fit theta over the whole nonlinear model. External backends
can define a more appropriate acceptance quantity without changing signal code.
Nonlinear scans currently require actual observations; background-only expected
limits are implemented for linear templates through each backend's own method.

No mutable global state, cache, implicit flux/density factor, silent normalization
or automatic convention conversion is introduced. Stored arrays/containers are
owned tuples; arrays passed into external providers/backends are copies. External
objects and callback closures must be pure, and their versions/configuration must
be recorded by the analyst. Returned counts and statistical metadata do not
replace source/calibration provenance.

## Example

```python
from wimp_legend.detector import DetectorResponse, ConstantQuenching, NoSmearing
from wimp_legend.backgrounds import ZeroBackground
from wimp_legend.statistics import PoissonCLsCounting
from wimp_legend.sensitivity import SensitivityAnalysis

analysis = SensitivityAnalysis(
    DetectorResponse(ConstantQuenching(1.), NoSmearing()),
    ZeroBackground(), PoissonCLsCounting(),
)
result = analysis.expected_limit(lambda E_nr_keV: 3., [0,.2,1,2], 5., (0,2))
# S=(3,12,15), B=(0,0,0), mu_up=2.302585092994046/30.
```

All values above are synthetic. Run `PYTHONPATH=src python
examples/generic_analysis.py` for quenching, Gaussian resolution, truth and
observed efficiencies, flat background, both expected-limit constructions, and
a disconnected nonlinear scan. Consult DETECTOR_RESPONSE.md, BACKGROUNDS.md and
STATISTICS.md for units, probability/exposure ownership and approximation limits.
