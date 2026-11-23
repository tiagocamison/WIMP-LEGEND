# Generic binned nuisance profiling

`ProfiledBinnedPoisson` is an opt-in statistics backend with a callback API.
It fits a nonnegative signal strength and any number of named nuisance
parameters. Existing `BinnedPoissonProfile` and `PoissonCLsCounting` calls
and conventions are unchanged. No physics or analysis-window selection is
performed: supply observations and expectations for exactly the bins to fit,
including sidebands or the whole configured window if desired.

```python
from wimp_legend.statistics import (
    GaussianPenalty, NuisanceParameter, ProfiledBinnedPoisson,
)

def expected_counts(mu, theta):
    scale = theta["scale"]
    return [2 * mu + 5 * scale, mu + 12 * scale]

model = ProfiledBinnedPoisson(
    expected_counts,
    [NuisanceParameter("scale", nominal=1, initial=.9, bounds=(0, None))],
    penalty=GaussianPenalty.independent({"scale": 1}, {"scale": .2}),
)
observed = [6, 11]
global_fit = model.fit(observed)
conditional_fit = model.fit(observed, mu=1)
observed_limit = model.observed_upper_limit(observed).mu_upper
median_asimov_limit = model.expected_upper_limit().mu_upper
```

Use `penalty=None` for an unconstrained fit, or provide a callable of the
read-only named nuisance mapping. The return value is in **−2 log likelihood
units**: an independent Gaussian adds `((theta-center)/sigma)**2`, with no
factor of one half. `GaussianPenalty(center, covariance)` accepts a symmetric
positive-definite covariance in the insertion order of the center mapping.
The Gaussian center and the nominal Asimov point are explicit, separate
inputs; normally use the same values. Expected limits use fractional
background-only counts evaluated at the nominal point. Observed fits
require integer counts.

Each fit returns the fitted values, penalized Poisson deviance (the requested
objective minus a data-only constant), evaluation count, optimizer status, and
a convergence flag. Equal lower and upper nuisance bounds fix a parameter.
The callbacks must be deterministic. Expected counts must be finite and
nonnegative; return `None` for points outside the model domain or use an
infinite penalty for excluded nuisance values. Put known domain constraints
in parameter bounds when possible. Supply a feasible configured start or
extra `starts=[(mu, {"name": value, ...}), ...]` for difficult fits.
The bounded Powell searches are local; repeat from physically plausible
seeds for multimodal models. `locally_flat` flags nearly constant coordinate
probes, not all correlated ridges or weak directions. A fitted nuisance
value on a ridge need not be unique. Failed fits are reported and unsuitable
for limits; a missing upper crossing is an error, not a finite exclusion.

Limits use the existing one-sided asymptotic convention:
`q_mu = 0` below the best-fit signal strength and
`q_mu = z_CL**2` at the upper crossing. This is **not CLs**.
The root finder assumes a single upper crossing, and the asymptotic
calibration is not validated for sparse counts, boundaries, or weak
identifiability. No coverage claim is made.

## Tester handoff

Run `python -m pytest -q -W error tests/test_statistics.py tests/test_statistics_nuisance.py`.
Check that the analytic two-bin fit recovers the conditional and unconditional
solutions, the correlated Gaussian gives the stated Mahalanobis penalty, zero
nuisances reproduce the legacy observed/Asimov limits, and a nuisance that
fully absorbs the signal cannot produce a finite limit. For future physical
callbacks, independently validate fit convergence from multiple starts and
the profile shape across the tested signal-strength range before using a
reported limit. Coverage and physical nuisance propagation are out of scope.
