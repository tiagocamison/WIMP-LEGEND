# Count-only statistical inference

Statistics consumes bin counts, not energy spectra, detector objects, WIMP
parameters or cross sections. Expected S/B vectors must be nonempty, finite,
nonnegative and have identical one-dimensional shapes. Actual observed n must
match and contain nonnegative integers (integer-valued floats are accepted).
There is no background subtraction. Actual observed counts/totals >=2^53 are
rejected because this float implementation cannot preserve integer precision.
No physical likelihood, data, nuisance prior or collaboration method is supplied.

## Public API (wimp_legend.statistics)

| API | Meaning |
| --- | --- |
| `StatisticsBackend.observed_upper_limit(signal,background,observed,*,confidence_level=.90)` | supplied templates and actual observations -> `LimitResult` |
| `StatisticsBackend.expected_upper_limit(signal,background,*,confidence_level=.90)` | background-only expected limit, construction defined by backend |
| `StatisticsBackend.acceptance_margin(signal,background,observed,*,confidence_level=.90)` | test absolute supplied signal at mu=1; nonnegative means accepted |
| `StatisticsBackend.acceptance_metadata()` | immutable `(method ID, construction, approximation)` for explicit parameter scans |
| `PoissonCLsCounting(root_policy=RootPolicy())` | exact known-background counting CLs; bins collapse to totals |
| `PoissonCLsCounting.log_cls(mu,signal,background,observed)` | natural log of the counting CLs ratio |
| `BinnedPoissonProfile(root_policy=RootPolicy())` | known-background independent-bin likelihood, asymptotic one-sided limits |
| `BinnedPoissonProfile.log_likelihood(mu,signal,background,observed)` | actual integer-count log likelihood |
| `BinnedPoissonProfile.fit(signal,background,observed)` | constrained maximum likelihood -> `FitResult` |
| `BinnedPoissonProfile.q_mu(mu,signal,background,observed)` | one-sided likelihood-ratio statistic |
| `RootPolicy(xtol=1e-10,rtol=1e-10,maxiter=200,max_bracket_steps=128)` | finite bracket growth followed by checked Brent root |
| `RootPolicy.crossing(function,lower=0,initial=1)` | increasing sign crossing; lower must be accepted/nonpositive; initial is a positive step |
| `LimitResult` | immutable limit and method metadata; fields below |
| `FitResult` | immutable mu_hat, log_likelihood, converged, status |

`LimitResult` records `mu_upper`, `method`, `construction`, `confidence_level`,
`kind`, `approximation`, `status`, `converged`, owned `reference_counts`, `notes`,
and `strength_convention`. The common strength convention is **lambda_i =
mu*S_i+B_i**, with S supplied at mu=1. The .90 default is visible in every
result. Infinity is reserved for `status='no_sensitivity'` when S is identically
zero. Bracket, precision and convergence failures raise errors, rather than
returning an apparently converged limit. Inputs are copied; results own tuples.
A collaboration backend can implement this protocol without subclassing either
fallback, and must return a LimitResult identifying its own construction.

## Exact Poisson CLs counting

Let S=sum S_i, B=sum B_i, n=sum n_i. The construction is

\[
\lambda=B+\mu S,\qquad
{\rm CL_s}(\mu)=\frac{P(N\le n\mid B+\mu S)}{P(N\le n\mid B)},\qquad
{\rm CL_s}(\mu_{\rm up})=1-CL.
\]

Known B is fixed and nonnegative. Shape information is deliberately discarded.
The result is an exact discrete Poisson CLs construction, with a numerical
root; it is not a classical central interval or a Gaussian S/sqrt(B) estimate.
Exact construction does not mean arbitrary-precision arithmetic or unqualified
coverage after including unmodeled background uncertainty.

For n=0 the background cancels analytically: CLs=exp(-mu*S), hence
`mu_up=-log(1-CL)/S`; S=1, CL=.90 gives **2.302585092994046**. This also works
with nonzero known B without subtracting large log-CDF exponents. More general
CDFs use scipy.stats.poisson.logcdf with a positive descending-series fallback
when it underflows. For B>n the common -B exponent is canceled analytically:

\[
\log CL_s=-s+n\log(1+s/B)+\log T_n(B+s)-\log T_n(B),\quad s=\mu S,
\]

where `T_n(lambda)=sum_{k=0}^n (n!/k!)*lambda^(k-n)`; the series is accumulated
from its largest term until further positive terms no longer change the sum.
This prevents erased CDF ratios in the tested large-background regime. Extremely
large counts near their mean can make the series expensive. Nonrepresentable
mean increments, strength results or totals are rejected. The root is solved in
signal **events**, then divided by S to reduce template-scale dependence.

The expected median is the distribution of limits under **N~Poisson(B)**.
Because the counting limit increases with n, take the **smallest integer n whose
Poisson CDF is >=.5**, and calculate that experiment's exact CLs limit. The
quantile bracket is checked. This is not n=B or a rounded Asimov observation.
The resulting kind is `expected_median`; reference_counts contains that one
collapsed count. Means >=2^52 are rejected for this quantile calculation.
Expected bands are not implemented in this milestone.

The scan acceptance margin is `log CLs(1)-log(1-CL)`.

## Independent-bin likelihood and asymptotic shape limits

\[
\lambda_i(\mu)=\mu S_i+B_i,\qquad
\log L(\mu)=\sum_i[n_i\log\lambda_i-\lambda_i-\log(n_i!)].
\]

For n=0, lambda=0 the term is zero. For n>0, lambda=0 the likelihood is zero
(log L=-infinity). If a positive observed bin has S=B=0, no strength has a
nonzero likelihood and fitting/limits raise ValueError.

Fit **mu_hat>=0** using the monotone likelihood score and a checked bracket.
A nonpositive score at zero fixes the fit to zero. For testing mu,

\[
q_\mu=\begin{cases}0,&\hat\mu>\mu,\\
2[\log L(\hat\mu)-\log L(\mu)],&\hat\mu\le\mu.\end{cases}
\]

The ratio is evaluated directly using log1p of mean increments, avoiding
subtraction of large absolute log likelihoods. Negative ratios from lost
precision are rejected, not clipped. The limit root above mu_hat satisfies

\[
q_{\mu_{\rm up}}=[\Phi^{-1}(CL)]^2.
\]

**CL must exceed .5** for this one-sided calibration. At .90,
z=Phi^-1(.90)=1.28155, and z^2 approximately 1.64237; this is not the
2.70554 two-sided chi-square threshold. Here Phi is the standard normal CDF.
This backend is **asymptotic only**. Results explicitly say
`status='asymptotic_unvalidated'` and note that **sparse-count quantitative
coverage is unsupported/not validated**. No threshold silently switches it to
another construction. High counts alone do not certify coverage for an arbitrary
analysis. Nuisance profiling, exact low-count shape calibration, Monte-Carlo CLs,
and uncertainty-aware collaboration likelihoods are future/external work.

Background-only expected shape sensitivity uses the **Asimov dataset n_i=B_i**,
including fractional values, only through `expected_upper_limit`. The kind is
`expected_asimov`; the actual-observation methods still reject fractional data.
This is an approximation, not a discrete median experiment. Counting and shape
results differ both in information retained and interval construction.

The scan acceptance margin is `z_CL^2-q_mu(1)`. For an arbitrary nonlinear
signal model, its alternative fit is along each tested point's fixed signal
shape; it is not a profile fit across that nonlinear parameter. A collaborator
must supply a different backend if that is the required test.

## Numerical scope and examples

RootPolicy grows a bracket by powers of two, verifies signs, calls Brent with
explicit tolerances/iteration budget, and checks convergence. No silent fallback
is applied after failure. Absolute tolerance is in the solver coordinate:
signal events for counting, mu for shape. Accuracy relative to a tiny root may
require a smaller explicit xtol. Very extreme float scales and cancellation can
still be unsupported; errors must be investigated. No arbitrary physical cutoff
or Python warning suppression is introduced.

```python
from wimp_legend.statistics import PoissonCLsCounting, BinnedPoissonProfile
exact = PoissonCLsCounting().observed_upper_limit([1.], [0.], [0])
median = PoissonCLsCounting().expected_upper_limit([1.], [3.2])
shape = BinnedPoissonProfile().expected_upper_limit([100.,0.], [1000.,9000.])
```

These are synthetic known-background templates. Retain the complete returned
metadata and supplied S/B/observations with any reported result.
