# WIMP-LEGEND independent detector-remediation re-audit

## 1. Verdict and audited source

**PASS WITH CONDITIONS for the supplied patched source tree. GA-AUDIT-001 and GA-AUDIT-002 are CLOSED.** Acceptance is bounded by the previously disclosed statistical and floating-point limits. This is not certification of a new remote commit: the fetched `refactor` tip still contained the failed-audit source, so the supplied patch was applied to the exact stated base for this audit.

| Source identity | Value |
| --- | --- |
| Repository | https://github.com/tiagocamison/WIMP-LEGEND |
| Exact base / fetched refactor tip | `166b1855ce75600bafb7a39947bd0887c31d8617` |
| Supplied patch | `CODE_AUDIT_003_REMEDIATION.diff` |
| Patch SHA-256 | `03a651bd564230d2a75e1f0747e296e9b5f9f5510f80c907e588f4019e8e73ed` |
| Exact resulting Git tree | `2b9447885e4b8c9d7cb09e42c2b41fc0b8d7efe5` |
| Audit date | 2026-10-19 |
| Python | 3.12.14, Clang 22.1.3 |
| NumPy / SciPy / pytest / mpmath | 2.3.5 / 1.17.0 / 9.1.1 / 1.4.1 |
| Installation | `python -m venv --system-site-packages audit4-venv`; editable package plus test dependencies and mpmath in that environment |

The base checkout was clean before applying the supplied patch. `git apply --check` succeeded. The patch was applied with `git apply --index`, and `git write-tree` recorded the immutable candidate identity before repair-source inspection. No auditor repair or change to production, inherited tests, or tolerances was made. At the end, the working tree matched the staged candidate: no unstaged changes, exactly the five supplied staged paths, and the same tree hash. A future commit is covered only if its audited files match this tree; the unpatched fetched branch is not approved by this verdict.

The supplied resolution record is byte-identical to the newly added `docs/CODE_AUDIT_003_RESOLUTION.md`. The patch changes only that record, detector documentation, `detector/response.py`, `detector/resolution.py`, and a new 20-case regression file. Rate, interactions, targets, halo, backgrounds, statistics, sensitivity, dependencies, and existing tests are unchanged.

## 2. Scope and independent evidence

The task is a remediation recheck of the previous generic-analysis audit, not a physical-data milestone. The audit:

1. recovered all prior reports and reconstructed their four test files without editing;
2. verified their published SHA-256 identities;
3. froze a reference/test plan before reading the repair implementation;
4. inspected every supplied change and its claimed mechanism;
5. ran the repository suite and all unchanged prior audit suites with warnings treated as errors;
6. added independent mathematical, property, range, and full-pipeline checks.

The new expected values use 80-digit mpmath erf/erfc expressions, elementary integrals, independent Lindhard inversion, and the previously written independent Poisson/likelihood oracles. No production detector/background/statistics method generates the expected values. The coder resolution and new tests were treated as claims and additional regression evidence, not as independent certification.

No physical Ge nuclear response, isotope composition, experimental calibration, LEGEND background, or measured exposure was introduced. No remote push, merge, or history rewrite was performed.

## 3. GA-AUDIT-001 recheck: declared efficiency boundaries

**CLOSED.** The repair now includes declared observed-efficiency boundaries, as well as observed bin edges, when constructing known truth-space crossings for `ConstantQuenching` and `LindhardQuenching`. Constant Gaussian widths also insert the same ±1,2,4,8 sigma neighborhoods around those boundaries. These are quadrature hints, not probability truncation or renormalization.

The original three unchanged counterexamples now pass:

- unit flat truth on `[0,1]`, observed acceptance only on `(0.4,0.41)`;
- `NoSmearing()`, `GaussianResolution(0)`, and `GaussianResolution(1e-5)`;
- accepted count 0.01 rather than the previous false zero.

New checks use nonunit q=0.07,0.4,1, disjoint acceptance windows, zero/narrow/broader Gaussian widths, and independent observed-space integrals. They also test deterministic Lindhard crossings for `(Z,k)=(8,0.1),(32,0.18),(80,0.3)`. Their expected truth crossings are computed from a separately written high-precision Lindhard formula, not from the production quenching method. A polynomial truth spectrum confirms that crossing placement preserves the integral rather than merely recovering a width.

Boundary hints are copied, deduplicated, and filtered. Caller arrays remain unchanged. The documented condition remains: arbitrary custom quenchers, callable widths, or custom kernels need their own truth integration hints. The repair does not claim automatic discovery of unreported callback features.

## 4. GA-AUDIT-002 recheck: local Gaussian probability coordinates

**CLOSED within the documented floating-point range.** Efficiency-weighted observed bins are now split in observed energy at declared boundaries and, where relevant, the Gaussian mean. Each local segment uses a CDF coordinate on the left or survival coordinate on the right. The inverse no longer receives a rounded probability of one merely because a coarse bin also spans the distribution's center. Segment contributions use `math.fsum`.

The original tail-window references are recovered:

| Mean, sigma=0.01 | Independent accepted probability for `(0.4,0.41)` |
| ---: | ---: |
| 0.30 | `7.619661958202802e-24` |
| 0.32 | `6.219831985865741e-16` |
| 0.33 | `1.2791904478283997e-12` |

The unchanged tests verify that the covering bin and a coarse `[0,1]` bin give the same finite result. The full discontinuous-efficiency fold with explicit truth hints also returns 0.01 instead of raising at an infinite efficiency coordinate.

Additional checks cover left tails, right tails, segments crossing the mean, several widths, arbitrary compatible bin partitions, duplicate/outside hints, callable-constant versus constant-efficiency shortcuts, and finite callback coordinates. Smooth affine efficiencies are checked against

\[
(\alpha+\beta m)p+\beta\sigma[\phi(z_a)-\phi(z_b)].
\]

Input/output scalar, empty, matrix and array shapes are checked. Invalid callback outputs, callback exceptions, and IntegrationWarning continue to fail explicitly. There is no clipping or window renormalization.

The helper includes a log-probability inverse when local interpolation underflows despite a positive segment mass. Source inspection verifies its mathematical orientation. This is not a claim that every branch is naturally reachable with the tested SciPy CDF implementation: `ndtr` can itself erase sufficiently remote tails before this fallback is reached. That retained boundary is demonstrated in Section 7.

The new `ndtri_exp` import is present in the [official SciPy 1.9.0 API documentation](https://docs.scipy.org/doc/scipy-1.9.0/reference/generated/scipy.special.ndtri_exp.html), so it does not by itself conflict with the project's `scipy>=1.9` declaration. Runtime execution here used SciPy 1.17.0; no cross-version numerical certification is claimed.

## 5. Regression results and source integrity

All test commands use `-q -W error`.

| Suite | Result |
| --- | --- |
| Patched repository | **675 passed in 67.89 s** |
| Unchanged pre-rate, rate, and generic audit files | **515 passed in 146.83 s** |
| New remediation extension | **50 passed, 1 failed in 224.31 s**; the sole failure is the disclosed remote-tail stress probe |

The repository count is 655 inherited cases +20 coder regression cases. The inherited external count is 154 pre-rate +237 rate +124 generic =515. The old seven detector falsifiers are retained unchanged inside those 124 generic cases; none is skipped or xfailed. The rate-product false-zero test remains part of the 237 rate cases.

The four recovered test hashes match the published records:

| Script | SHA-256 |
| --- | --- |
| `test_independent.py` | `1f1827d261b49a0928db7eec5c33c0bc5d518d976d107f245fd4de6639ac23e6` |
| `test_falsifications.py` | `99438a0db818631638f99d80b7282c8f6ed8d56505443abf4057b4e0005de6cb` |
| `test_rate_audit.py` | `a90623aca9133d9293181cad00c4db5dbbb5fbc0f72cc72793612c172ba24aae` |
| `test_generic_audit.py` | `1705ebe6e24403a54573d3dd6b0ba8724c681bbc29a1e684a93e8bdc1fc05aec` |

Only the supplied candidate patch modifies the source checkout. `git diff` against the index is empty; the staged candidate remains exactly the identified tree. The frozen NREFT convention remains unchanged. The new independent test file and its strict remote-tail assertion are included verbatim below.

## 6. Independent discontinuous-efficiency pipeline

The new pipeline differs from both the original smooth audit benchmark and the coder's repair benchmark:

- truth density 3 events/(kg day keV_nr), truth range `[0,3]` keV_nr;
- constant q=0.4, Gaussian sigma=0.02 keV;
- exposure 17 kg day;
- observed edges `(0,0.41,0.9,1.2)` keV;
- observed efficiency 0.4 on `(0.22,0.28)`, 0.8 on `(0.7,0.83)`, zero elsewhere;
- differential background `b(y)=0.3+0.2y` events/(kg day keV_obs);
- actual observations `(1,2,0)`, CL=0.90.

Swap the truth and observed integrations. For flat truth density r and finite truth interval `[0,L]`, the independent observed differential density is

\[
g(y)=\frac rq\left[\Phi(y/\sigma)-\Phi((y-qL)/\sigma)\right].
\]

Integrating this over the actual acceptance windows gives S; a quadratic primitive gives B. The independent Poisson finite-sum and likelihood oracles then give the four limits, without calling production detector/background/statistics methods for expected values.

| Bin | Independent S | Independent B | Observed n |
| --- | ---: | ---: | ---: |
| [0,0.41) | 3.060000000000002 | 2.3767699999999996 | 1 |
| [0.41,0.9) | 13.260000000000005 | 3.590230000000001 | 2 |
| [0.9,1.2] | 0 | 2.6009999999999995 | 0 |

The total background is 8.568 and the independent discrete median observation is 8. The zero signal in the third bin is an actual acceptance choice, not a missed truth feature.

| Construction | Independent mu upper |
| --- | ---: |
| Exact counting, observed | 0.1923154443637651 |
| Exact counting, discrete median expected | 0.35616277791362394 |
| Asymptotic shape, observed | 0.09372825543987053 |
| Asymptotic shape, expected Asimov | 0.2137150370645395 |

The pipeline is executed with default integration/root settings and with `epsrel=1e-10`, `xtol=rtol=1e-13`. Intermediate S/B and all limits are checked separately. The execution evidence in Appendix C records both sets of actual limits. The previous smooth pipeline is also retained inside the unchanged generic suite.

Maximum new-pipeline signal discrepancy is approximately `5.33e-15` events (about `4.02e-16` relative in the nonzero second bin); background counts agree to the reported digits. Default-versus-tight limit differences are at most about `4.97e-12` in mu. All four independent limit assertions and both tolerance-convergence tests pass.

The implementation's own 20 added cases supply useful independent normal/Poisson/likelihood regressions, including a different discontinuous pipeline. Some coder quenching-boundary tests construct their observed boundaries with the production quenching model; those are consistency properties rather than independent tests of the model formula. The new Lindhard checks above remove that shared-oracle dependence.

## 7. Remaining conditions and numerical stress result

### GA-AUDIT-003 — retained NOTE: asymptotic shape coverage

The shape backend remains explicitly asymptotic, with `asymptotic_unvalidated` metadata and no claim of calibrated sparse-count coverage. Statistics code is byte-identical to the prior audit. The prior sparse-toy study need not be rerun to validate a detector-only repair. Its demonstrated undercoverage remains an acceptance condition for scientific use. Exact counting and asymptotic shape remain different constructions.

### GA-AUDIT-004 — retained NOTE: root coordinate and template scaling

Shape root absolute tolerance remains in mu, while counting roots are solved in signal events. Very small strength roots need appropriate template normalization and an explicit tolerance/convergence check. The unchanged tighter-tolerance rescaling assertion is rerun. No new statistical calibration or solver fix is claimed.

### RM-AUDIT-001 — NOTE: documented remote Gaussian tails can still be erased

- **Location:** `src/wimp_legend/detector/resolution.py`, `GaussianResolution.bin_probabilities` and `_weighted_bin_probability`, CDF/survival mass evaluation before interpolation.
- **Expected under a strict all-representable-values criterion:** a finite interval with nonzero representable probability remains nonzero, or produces an explicit unsupported-range error.
- **Observed:** in the new piecewise-partition stress case with mean=0.15, sigma=0.009, the interval `[0.49,0.83]` begins at about 37.78 sigma. Weighted probability 0.4 times its independent erfc mass is approximately `5.26582285e-313`; the implementation returns zero. This is a subnormal tail, not an ordinary acceptance probability.
- **Reproducer:** `test_local_piecewise_partition[0.009-0.15]` in Appendix A.
- **Interpretation:** the assertion deliberately applies the stricter all-representable-values property and remains failed in the reported stress results. The resolution record and detector documentation explicitly retain extremely remote-tail/intermediate-underflow limits. The previous generic audit already documented the same class at 38 sigma. This is therefore a reconfirmed scope boundary, not a new regression or an unclosed ordinary-scale MAJOR issue.
- **Scientific consequence:** do not rely on extreme far-tail migration or arbitrary rescaling of its probability into large event counts without a separate supported numerical analysis. The supplied repair does not offer arbitrary-precision tail masses.
- **Recommendation:** retain this qualification and the probe. If such tails become part of the intended scientific range, compute local segment mass in stable log/tail arithmetic and enforce a quantitative supported range or explicit failure. The log inverse for interpolation alone cannot restore a mass already set to zero by the CDF.

**Severity accounting:** no unresolved BLOCKER, MAJOR, or MINOR is identified within the declared remediation scope. Three NOTE-level conditions remain: two inherited statistical conditions and the reconfirmed remote-tail bound. The two old MAJOR findings are closed. The strict stress failure is disclosed and is not relabeled as a passing assertion, skipped, or xfailed.

## 8. Subsystem and milestone verdicts

| Subsystem | Verdict | Basis / condition |
| --- | --- | --- |
| Exact source/evidence identity | PASS for the supplied patch tree | Exact base, patch hash, and resulting tree; no new remote commit supplied |
| GA-AUDIT-001 closure | CLOSED / PASS | Declared efficiency boundaries reach truth quadrature for built-in monotone models |
| GA-AUDIT-002 closure | CLOSED / PASS WITH CONDITIONS | Local tail coordinates fix original coarse-bin zero/infinity; remote-tail mass limit retained |
| Quenching and units | PASS | Unchanged implementation; independent formulas and inversion checks pass |
| Ordinary detector migration/efficiency | PASS | Original seven counterexamples and independent partition/affine/pipeline checks |
| Backgrounds | PASS inherited and rechecked | Unchanged source; new full-pipeline event counts checked independently |
| Exact counting statistics | PASS inherited and rechecked | Unchanged source; independent four-limit pipeline |
| Asymptotic shape statistics | PASS WITH CONDITIONS | Sparse coverage and solver-coordinate qualification unchanged |
| Generic orchestration/external ports | PASS inherited | Source unchanged; prior external-port and nonlinear tests rerun |
| Full discontinuous pipeline | PASS | Default/tighter integration and roots, intermediate counts, both backends |
| Architecture and data containment | PASS | Detector-only repair, no new reverse model dependency or physical data |

**Overall: PASS WITH CONDITIONS for tree `2b9447885e4b8c9d7cb09e42c2b41fc0b8d7efe5`.** The supplied generic analysis foundation is accepted to proceed to a **separately sourced and validated physical Ge data milestone**, within the documented numerical/statistical scope. This does not certify any nuclear coefficients, calibration, composition, background model, uncertainty treatment, or experimental sensitivity. It also does not approve the unpatched remote branch.

Before importing physical data, preserve a commit containing exactly this accepted repair tree and its base/patch identity; keep the unchanged audit chain and new probes reproducible. Physical nuclear conventions, provenance, source momentum mapping/ranges, missing-versus-zero channels, target counts and detector calibration still need their own validation. Use a statistically appropriate backend and converged numerical policies for the intended data regime. No physical data import was performed in this audit.

## Appendix A. New independent remediation tests

Save verbatim as `audit4-work/test_remediation.py`. This imports only the mathematical oracles `cls_root`, `median`, and `shape_ref` from the unchanged independent `test_generic_audit.py`, not repository tests. The strict remote-tail failing assertion is deliberately retained.

```python
"""Independent remediation extension: expected values from mathematical kernels."""
import json,math
from pathlib import Path
from dataclasses import replace
import mpmath as mp
import numpy as np
import pytest
from scipy.integrate import quad,IntegrationWarning
from scipy.optimize import brentq
from wimp_legend.detector import (ConstantQuenching,LindhardQuenching,NoSmearing,GaussianResolution,CallableEfficiency,ConstantEfficiency,DetectorResponse,IntegrationPolicy)
from wimp_legend.backgrounds import CallableBackground
from wimp_legend.statistics import PoissonCLsCounting,BinnedPoissonProfile,RootPolicy
from wimp_legend.sensitivity import SensitivityAnalysis
# Imported independent audit ORACLES, not repository tests or production expected generators.
from test_generic_audit import cls_root,median,shape_ref
mp.mp.dps=80

def interval(a,b):
 a,b=mp.mpf(float(a)),mp.mpf(float(b));s=mp.sqrt(2)
 if a>=0:return float((mp.erfc(a/s)-mp.erfc(b/s))/2)
 if b<=0:return float((mp.erfc(-b/s)-mp.erfc(-a/s))/2)
 return float((mp.erf(b/s)-mp.erf(a/s))/2)

def lindhard(E,Z,k):
 E=mp.mpf(float(E));x=mp.mpf('11.5')*E*mp.mpf(Z)**(-mp.mpf(7)/3)
 g=3*x**mp.mpf('.15')+mp.mpf('.7')*x**mp.mpf('.6')+x
 return float(E*mp.mpf(float(k))*g/(1+mp.mpf(float(k))*g))

@pytest.mark.parametrize('mean',[0.,.15,.43,.8,1.2])
@pytest.mark.parametrize('sigma',[.009,.17,2.])
def test_local_piecewise_partition(mean,sigma):
 boundaries=np.array([-.3,.21,.22,.49,.83,1.4]);heights=[.1,.7,0.,.4,.8]
 def eff(y):
  assert np.isfinite(y) and boundaries[0]-1e-13<=y<=boundaries[-1]+1e-13
  return heights[min(max(np.searchsorted(boundaries,y,side='right')-1,0),4)]
 ref=np.array([h*interval((a-mean)/sigma,(b-mean)/sigma) for h,a,b in zip(heights,boundaries[:-1],boundaries[1:])])
 g=GaussianResolution(sigma)
 hints=np.array([.83,.21,.22,.49,.21,-10.,10.]);saved=hints.copy()
 coarse=g.bin_probabilities(mean,[boundaries[0],boundaries[-1]],observed_efficiency=CallableEfficiency(eff),observed_breakpoints=hints)[0]
 fine=g.bin_probabilities(mean,boundaries,observed_efficiency=CallableEfficiency(eff))
 np.testing.assert_allclose(fine,ref,rtol=8e-12,atol=0)
 assert coarse==pytest.approx(math.fsum(ref),rel=8e-12,abs=0)
 np.testing.assert_array_equal(hints,saved)

@pytest.mark.parametrize('mean,bins',[(0.,[-8,-7.9]),(0.,[7.9,8]),(.4,[-.2,1.]),(2.,[1.99,2,2.01]),(0.,[37.5,37.6])])
def test_callable_constant_matches_shortcut(mean,bins):
 model=GaussianResolution(1)
 expected=np.array([.4*interval(a-mean,b-mean) for a,b in zip(bins[:-1],bins[1:])])
 fn=model.bin_probabilities(mean,bins,observed_efficiency=CallableEfficiency(lambda y:.4))
 direct=model.bin_probabilities(mean,bins,observed_efficiency=ConstantEfficiency(.4))
 np.testing.assert_allclose(fn,expected,rtol=3e-11,atol=0)
 np.testing.assert_allclose(direct,expected,rtol=3e-11,atol=0)

@pytest.mark.parametrize('mean',[0.,.35,1.])
def test_affine_kernel_accuracy(mean):
 sig=.2;bins=[-.2,.2,.7,1.4]
 phi=lambda z:math.exp(-z*z/2)/math.sqrt(2*math.pi)
 ref=[]
 for a,b in zip(bins[:-1],bins[1:]):
  za,zb=(a-mean)/sig,(b-mean)/sig
  ref.append((.5+.1*mean)*interval(za,zb)+.1*sig*(phi(za)-phi(zb)))
 got=GaussianResolution(sig).bin_probabilities(mean,bins,observed_efficiency=CallableEfficiency(lambda y:.5+.1*y),policy=IntegrationPolicy(epsrel=1e-11))
 np.testing.assert_allclose(got,ref,rtol=3e-11,atol=2e-15)

@pytest.mark.parametrize('Z,k',[(8,.1),(32,.18),(80,.3)])
@pytest.mark.parametrize('width',[None,0.])
def test_lindhard_crossings_from_independent_formula(Z,k,width):
 a,b=.1,.103
 l=brentq(lambda E:lindhard(E,Z,k)-a,0,10,xtol=1e-14)
 u=brentq(lambda E:lindhard(E,Z,k)-b,0,10,xtol=1e-14)
 primitive=lambda E:E+.2*E*E
 ref=primitive(u)-primitive(l)
 d=DetectorResponse(LindhardQuenching(Z,k),NoSmearing() if width is None else GaussianResolution(width),observed_efficiency=CallableEfficiency(lambda y:1. if a<y<b else 0.),observed_breakpoints=[a,b])
 assert d.expected_counts(lambda E:1+.4*E,[0,2],1,(0,10))[0]==pytest.approx(ref,rel=2e-11)

@pytest.mark.parametrize('q',[.07,.4,1.])
@pytest.mark.parametrize('sigma',[0.,1e-5,.03])
def test_nonunit_quenching_disjoint_acceptance(q,sigma):
 # Acceptance edges expressed in observed energy from independently chosen truth
 # fractions; Gaussian reference swaps integrals instead of using detector kernels.
 L=3.;a,b=q*.7,q*.76;d,e=q*1.6,q*1.9
 eps=lambda y:.4 if a<y<b else (.8 if d<y<e else 0.)
 model=DetectorResponse(ConstantQuenching(q),GaussianResolution(sigma),observed_efficiency=CallableEfficiency(eps),observed_breakpoints=[a,b,d,e])
 if sigma==0:ref=(.4*(b-a)+.8*(e-d))/q
 else:
  def density(y):return interval((y-q*L)/sigma,y/sigma)/q
  ref=.4*quad(density,a,b,epsabs=1e-12,epsrel=1e-11)[0]+.8*quad(density,d,e,epsabs=1e-12,epsrel=1e-11)[0]
 actual=model.expected_counts(lambda E:1,[0,q*L],1,(0,L))[0]
 assert actual==pytest.approx(ref,rel=2e-9,abs=0)

@pytest.mark.parametrize('shape',[(),(0,),(2,0),(2,2)])
def test_shapes_and_finite_interval_coordinates(shape):
 E=np.full(shape,.31)
 model=GaussianResolution(.02)
 before=E.copy()
 def eff(y):
  assert np.isfinite(y) and .2-1e-14<=y<=.7+1e-14
  return .3
 out=model.bin_probabilities(E,[.2,.4,.7],observed_efficiency=CallableEfficiency(eff))
 assert np.shape(out)==shape+(2,)
 np.testing.assert_array_equal(E,before)

@pytest.mark.parametrize('value',[-.1,1.1,np.nan,np.inf,[.2],True])
def test_invalid_weighted_callbacks_remain_errors(value):
 with pytest.raises((ValueError,TypeError)):
  GaussianResolution(.1).bin_probabilities(.4,[0,1],observed_efficiency=lambda y:value,observed_breakpoints=[.3,.7])

def test_warning_and_callback_exceptions_propagate(monkeypatch):
 import warnings,wimp_legend._analysis as shared
 def broken(y):raise RuntimeError('external failure')
 with pytest.raises(RuntimeError):GaussianResolution(.1).bin_probabilities(.4,[0,1],observed_efficiency=broken)
 def fail(*a,**k):warnings.warn('independent audit injected warning',IntegrationWarning);return 0.,0.
 monkeypatch.setattr(shared,'quad',fail)
 with pytest.raises(IntegrationWarning):GaussianResolution(.1).bin_probabilities(.4,[0,1],observed_efficiency=lambda y:.2)

@pytest.mark.parametrize('backend_kind',['exact','shape'])
def test_discontinuous_pipeline_independent_four_limits(backend_kind):
 edges=np.array([0.,.41,.9,1.2]);q=.4;sig=.02;X=17.;L=3.
 windows=[(.22,.28,.4),(.7,.83,.8)]
 def acceptance(y):
  return sum(w for a,b,w in windows if a<y<b)
 def density(y):return 3/q*interval((y-q*L)/sig,y/sig)
 S=[]
 for lo,hi in zip(edges[:-1],edges[1:]):
  s=0.
  for a,b,w in windows:
   lower,upper=max(lo,a),min(hi,b)
   if lower<upper:s+=w*quad(density,lower,upper,epsabs=1e-12,epsrel=1e-11)[0]
  S.append(X*s)
 S=np.array(S);B=X*np.diff(.3*edges+.1*edges*edges);n=np.array([1,2,0])
 refs=[cls_root(sum(S),sum(B),int(sum(n))),cls_root(sum(S),sum(B),median(sum(B)))] if backend_kind=='exact' else [shape_ref(S,B,n)[1],shape_ref(S,B,B)[1]]
 results=[]
 for tight in [False,True]:
  policy=IntegrationPolicy(epsrel=1e-10 if tight else 1e-8)
  roots=RootPolicy(xtol=1e-13 if tight else 1e-10,rtol=1e-13 if tight else 1e-10)
  backend=PoissonCLsCounting(roots) if backend_kind=='exact' else BinnedPoissonProfile(roots)
  detector=DetectorResponse(ConstantQuenching(q),GaussianResolution(sig),observed_efficiency=CallableEfficiency(acceptance),observed_breakpoints=[.22,.28,.7,.83],policy=policy)
  analysis=SensitivityAnalysis(detector,CallableBackground(lambda y:.3+.2*y),backend)
  observed=analysis.observed_limit(lambda E:3.,edges,X,(0,L),n)
  expected=analysis.expected_limit(lambda E:3.,edges,X,(0,L))
  results.append([observed.limit.mu_upper,expected.limit.mu_upper])
  for actual,ref in zip([observed,expected],refs):
   np.testing.assert_allclose(actual.signal_counts,S,rtol=2e-10,atol=1e-14)
   np.testing.assert_allclose(actual.background_counts,B,rtol=2e-14)
   assert actual.limit.mu_upper==pytest.approx(ref,rel=3e-9)
 np.testing.assert_allclose(results[0],results[1],rtol=3e-9)
 Path(__file__).with_name('pipeline_'+backend_kind+'.json').write_text(json.dumps(dict(S_reference=S.tolist(),B_reference=B.tolist(),observed=n.tolist(),reference_limits=refs,actual_limits=results,signal_actual=observed.signal_counts,background_actual=observed.background_counts),indent=2))
```

## Appendix B. Pre-inspection reference plan

```text
Remediation reference plan frozen before repair-source inspection.
Base 166b1855ce75600bafb7a39947bd0887c31d8617 plus supplied patch; record immutable tree before tests.
Inherited tests must remain byte-identical: 154 pre-rate,237rate,124generic. Main falsifiers: observed acceptance (.4,.41), unit truth0..1 gives .01; local weighted-tail bin probabilities equal erf/erfc interval irrespective of coarse/fine binning. Quenching q changes deterministic accepted truthwidth by1/q. Both constant and monotone Lindhard map declared efficiency edges to truth crossings. Finite Gaussian integration should retain tailmass, not truncate or renormalize.
Extend tests: left/right/crossing-mean segments; affine and discontinuous efficiencies; many/duplicate/outside boundaries; scalar/array/empty coordinates; nonunit constant q and multiple Lindhard models; partition additivity; constant-efficiency and callable-constant equivalence in central/tail bins; subnormal local probability interpolation; callback energies finite/in interval; epsrel convergence and unchanged callback errors/warnings.
Independent formulas: p(a,b)=.5*(erfc((a-m)/sqrt(2)/s)-erfc((b-m)/sqrt(2)/s)); use symmetry for left tail and mpmath80digits. Weighted affine (alpha+beta*m)*p+beta*sigma*(phi(za)-phi(zb)). Piecewise-constant efficiency sum(weights*p) over actual represented segment endpoints. Truth unit density, constant q, finite rangeL..U: observed density=[Phi((y-qL)/sigma)-Phi((y-qU)/sigma)]/q. Swap integration to observed y; analytic normal primitive F(z)=zPhi(z)+phi(z) integrates this exactly for acceptance windows. Polynomial truth tests independent nested integration or integration by parts. Matrix/kernel/stats/rate unchanged verified by diff.
Pipeline independently reuse old highprecision likelihood/CLs oracles only, own different detector fixture: truth density3, q=.4, sigma=.02, truth[0,3], exposure17, observed edges[0,.41,.9,1.2], disjoint observed acceptance (.22,.28) weight.4 and(.7,.83) weight.8; background b=.3+.2y, observed n=[1,2,0]. Derive S by swapped finite-truth Gaussian integral; B by primitive. Four limits independent inherited mathematical refs; no production expected generator. Compare default/tight quadrature and root policies.
No physical data or production repair by auditor. Conserve existing sparse-count asymptotic and root-coordinate conditions.
```

## Appendix C. Complete execution evidence and reproduction

Exact new-suite failure evidence:

```text
=================================== FAILURES ===================================
__________________ test_local_piecewise_partition[0.009-0.15] __________________

mean = 0.15, sigma = 0.009

    @pytest.mark.parametrize('mean',[0.,.15,.43,.8,1.2])
    @pytest.mark.parametrize('sigma',[.009,.17,2.])
    def test_local_piecewise_partition(mean,sigma):
     boundaries=np.array([-.3,.21,.22,.49,.83,1.4]);heights=[.1,.7,0.,.4,.8]
     def eff(y):
      assert np.isfinite(y) and boundaries[0]-1e-13<=y<=boundaries[-1]+1e-13
      return heights[min(max(np.searchsorted(boundaries,y,side='right')-1,0),4)]
     ref=np.array([h*interval((a-mean)/sigma,(b-mean)/sigma) for h,a,b in zip(heights,boundaries[:-1],boundaries[1:])])
     g=GaussianResolution(sigma)
     hints=np.array([.83,.21,.22,.49,.21,-10.,10.]);saved=hints.copy()
     coarse=g.bin_probabilities(mean,[boundaries[0],boundaries[-1]],observed_efficiency=CallableEfficiency(eff),observed_breakpoints=hints)[0]
     fine=g.bin_probabilities(mean,boundaries,observed_efficiency=CallableEfficiency(eff))
>    np.testing.assert_allclose(fine,ref,rtol=8e-12,atol=0)
E    AssertionError: 
E    Not equal to tolerance rtol=8e-12, atol=0
E    
E    Mismatched elements: 1 / 5 (20%)
E    Max absolute difference among violations: 5.26582285e-313
E    Max relative difference among violations: 1.
E     ACTUAL: array([1.000000e-01, 9.156164e-12, 0.000000e+00, 0.000000e+00,
E           0.000000e+00])
E     DESIRED: array([1.000000e-001, 9.156164e-012, 0.000000e+000, 5.265823e-313,
E           0.000000e+000])

audit4-work/test_remediation.py:41: AssertionError
=========================== short test summary info ============================
FAILED audit4-work/test_remediation.py::test_local_piecewise_partition[0.009-0.15]
1 failed, 50 passed in 224.31s (0:03:44)
```

Pipeline exact evidence (actual limits are default, then tighter policies):

```json
{
  "S_reference": [
    3.060000000000002,
    13.260000000000005,
    0.0
  ],
  "B_reference": [
    2.3767699999999996,
    3.590230000000001,
    2.6009999999999995
  ],
  "observed": [
    1,
    2,
    0
  ],
  "reference_limits": [
    0.1923154443637651,
    0.35616277791362394
  ],
  "actual_limits": [
    [
      0.1923154443637675,
      0.35616277791447104
    ],
    [
      0.1923154443637675,
      0.35616277791362405
    ]
  ],
  "signal_actual": [
    3.060000000000002,
    13.26,
    0.0
  ],
  "background_actual": [
    2.3767699999999996,
    3.590230000000001,
    2.6009999999999995
  ]
}
```

Pipeline shape evidence (actual limits are default, then tighter policies):

```json
{
  "S_reference": [
    3.060000000000002,
    13.260000000000005,
    0.0
  ],
  "B_reference": [
    2.3767699999999996,
    3.590230000000001,
    2.6009999999999995
  ],
  "observed": [
    1,
    2,
    0
  ],
  "reference_limits": [
    0.09372825543987053,
    0.2137150370645395
  ],
  "actual_limits": [
    [
      0.09372825543980798,
      0.21371503705957437
    ],
    [
      0.09372825543980798,
      0.21371503706453954
    ]
  ],
  "signal_actual": [
    3.060000000000002,
    13.26,
    0.0
  ],
  "background_actual": [
    2.3767699999999996,
    3.590230000000001,
    2.6009999999999995
  ]
}
```

Additional source/evidence hashes:

| File | SHA-256 |
| --- | --- |
| `upload/CODE_AUDIT_003_REMEDIATION.diff` | `03a651bd564230d2a75e1f0747e296e9b5f9f5510f80c907e588f4019e8e73ed` |
| `upload/CODE_AUDIT_003_RESOLUTION.md` | `4093b63ffc0de8a78ab92a79318ee447e883344d0da5c6cecbd8fe2576be2e7a` |
| `audit4-work/test_remediation.py` | `99110bc9daf4165844d3439e3f46a6ce85e7acfb6b9086b912a13f4813dd1847` |
| `audit4-work/reference_plan.md` | `24a355ddb8a15ca5ee0bcb3b5e8f2800dcbc8eb01075deaf2fda2b1ff72f3382` |
| `audit-repo/src/wimp_legend/detector/resolution.py` | `0c4b0685efe85f1e7c6ff6b1f3f2715b2080cac8d244c8d3cd06542488d794e6` |
| `audit-repo/src/wimp_legend/detector/response.py` | `5264f71c0f1d55c12bff6ef768d3d524615fc88b7a787b550105fb62d5002cf6` |
| `audit-repo/docs/NREFT_CONVENTIONS.md` | `a2c57658add45b273ca082d8b01fbab7dd739f8625d6f04a438d619d4cf385cd` |

Reproduction from a workspace with the supplied patch in `upload/`:

```bash
git clone --branch refactor --single-branch https://github.com/tiagocamison/WIMP-LEGEND.git audit-repo
git -C audit-repo switch --detach 166b1855ce75600bafb7a39947bd0887c31d8617
git -C audit-repo apply --check ../upload/CODE_AUDIT_003_REMEDIATION.diff
git -C audit-repo apply --index ../upload/CODE_AUDIT_003_REMEDIATION.diff
git -C audit-repo write-tree
# Expected: 2b9447885e4b8c9d7cb09e42c2b41fc0b8d7efe5
python -m venv --system-site-packages audit4-venv
audit4-venv/bin/python -m pip install -e './audit-repo[test]' mpmath==1.4.1
audit4-venv/bin/python -m pytest -q -W error audit4-work/test_independent.py audit4-work/test_falsifications.py audit4-work/test_rate_audit.py audit4-work/test_generic_audit.py
audit4-venv/bin/python -m pytest -q -W error audit4-work/test_remediation.py
```

Run the repository suite from `audit-repo` with `../audit4-venv/bin/python -m pytest -q -W error`. Extract the first two unchanged test files from CODE_AUDIT_001_PRE_RATE.md Appendix A, the rate file from CODE_AUDIT_002_RATE_LAYER.md Appendix A, and the generic file from CODE_AUDIT_003_GENERIC_ANALYSIS.md Appendix A. Verify all four hashes before execution. The inherited architecture tests expect the source directory to be named `audit-repo` beside `audit4-work`.

The new suite's strict stress assertion exits nonzero because it probes beyond the documented ordinary tail range; this report does not instruct marking it xfail or altering its tolerance. The 675+515 inherited/regression cases all pass, while the 51 additional cases have the explicitly reported 50/1 result. Test timings are environmental.
