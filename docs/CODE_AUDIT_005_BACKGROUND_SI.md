# WIMP-LEGEND independent background-routing and SI-normalization audit

## 1. Verdict and exact source identity

**Background routing: PASS. SI normalization: PASS. Overall: PASS WITH CONDITIONS**, retaining the already documented numerical and statistical scope of the accepted foundation. No new implementation defect was identified in the two additions. This certificate applies to the verified delivered **source tree**, not independently authenticated commit metadata or physical experimental inputs.

| Source item | Identity |
| --- | --- |
| Repository | https://github.com/tiagocamison/WIMP-LEGEND |
| Verified base commit | `b48fc6d040ec30dc3c008a981773bba16bff300f` |
| Candidate commit declared by handoff | `9d2468e875792cd19210904e401d383d242d6cc5` |
| Verified candidate tree | **`2da8deb5cac7c9f07c73722f0288e1ced38a5c57`** |
| Delivered diff SHA-256 | `f2f237d0d94cbd2723d7d283bfa8982ca39733160d7247342d1deb06db6c72c3` |
| Frozen convention | `anand2014_prc89_065501_half_isospin_v1` |
| Audit date | 2026-10-19 |

The handoff explicitly permits applying the diff to the exact base and verifying `git write-tree` as an alternative to importing the Git bundle. The attachments contain the diff and two reports, but no bundle or SOURCE_IDENTITY.json. The fetched remote branch reaches the exact base. Applying the supplied diff with `git apply --index` gives the exact handoff tree above. Thus source-content identity is established; the declared candidate commit object, parent/committer metadata, and local branch identity are not independently recovered. No new commit was invented for certification.

The base includes commit `31d823f76de7635816db66238b55ec952428abfa`, whose tree `2b9447885e4b8c9d7cb09e42c2b41fc0b8d7efe5` exactly matches the detector repair accepted in CODE_AUDIT_004. The intervening base commit adds audit documentation. This preserves the accepted remediation rather than silently switching to the earlier failed detector source.

The candidate has twelve changed paths: routing/energy/ER/SI modules and exports, the detector point-line methods/domain declarations, two new test files, and three documentation files. All changes were inspected. The working tree matches the staged candidate, and its tree hash remains unchanged after testing. The auditor applied only the delivered patch and did not repair production code, alter repository tests, or change any inherited audit assertion or tolerance.

## 2. Independent method and scientific scope

The frozen NREFT convention and prior independent audits are the analytic and validation context. Before inspecting new implementation or coder regressions, a reference note recorded the background continuum/line equations, NR and ER mappings, absolute SI conversion, half-isospin transformation, coherent cross section, polynomial-speed-PDF rate, and independent synthetic pipeline configurations. The note is included below.

The new audit then constructed its own synthetic masses/charges, responses, coefficients, lines and background spectra. Expected values use mpmath at 90 decimal digits, elementary polynomial integrals, erfc Gaussian probabilities, and explicit Poisson or one-sided likelihood relations. No production response/rate/normalization helper generates the expected values. Production methods are used for values under test. Coder tests were read after writing the independent references and new test file.

No physical Ge coefficient, empirical isotope mass, decay/neutrino model, calibration, LEGEND background, or observed experimental dataset was added. The proton normalization mass in the tests is an explicit synthetic parameter, not a newly imported physical mass. Missing physical inputs are not defects of this generic milestone and are not certified by this result.

### Environment and commands

Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0, pytest 9.1.1, mpmath 1.4.1. The existing isolated `audit4-venv` was reused; its editable package path points to the current `audit-repo` candidate. The preceding audit's source checkout was preserved separately. No project dependency specification was changed.

Every suite below runs with warnings treated as errors:

```bash
# From the candidate checkout:
../audit4-venv/bin/python -m pytest -q -W error
# From its workspace parent:
audit4-venv/bin/python -m pytest -q -W error audit4-work/test_independent.py audit4-work/test_falsifications.py audit4-work/test_rate_audit.py audit4-work/test_generic_audit.py
audit4-venv/bin/python -m pytest -q -W error audit4-work/test_remediation.py
audit4-venv/bin/python -m pytest -q -W error audit5-work/test_background_si_audit.py
```

Final execution totals and evidence are recorded in Section 6 and Appendix C. No skip or xfail was introduced.

## 3. A — Background routing

### Independent equations and units

For a continuum in the declared truth variable E,

\[
B_j=X\int r(E)\epsilon_{\rm truth}(E)
\int_{a_j}^{b_j}K(y\mid E)\epsilon_{\rm obs}(y)\,dy\,dE.
\]

Here X is kg day, r is events/(kg day keV_nr) or events/(kg day keV_er), and the result is an expected event count. With deterministic constant NR quenching q, the observed density is `r(y/q)/q`; ER identity mapping uses q=1. Already observed inputs are integrated or binned directly, without this response integral.

For a monoenergetic truth line of integrated intensity I in events/(kg day),

\[
B_{j,\rm line}=XI\,\epsilon_{\rm truth}(E_l)
\int_{a_j}^{b_j}K(y\mid E_l)\epsilon_{\rm obs}(y)dy.
\]

**There is no `1/q` factor in the integrated line area.** The Jacobian of a transformed delta measure cancels when integrated. Line input is not a sampled differential spike to be found by truth quadrature.

For Gaussian mean m and sigma s, expected bin probability is

\[
p(a,b)=\Phi((b-m)/s)-\Phi((a-m)/s),
\]

evaluated independently through high-precision erfc to retain tails. Affine observed efficiency `alpha+beta*y` uses the exact weighted probability

\[
(\alpha+\beta m)p+\beta s[\phi(z_a)-\phi(z_b)].
\]

### Results

| Requirement | Independent evidence | Verdict |
| --- | --- | --- |
| NR truth quenching and continuum normalization | Polynomial `r(E)=2+E`, q=0.2 and 0.65, explicit truth range, nonuniform observed bins and X=0,1,11; exact primitive with transformed limits | PASS |
| ER bypasses NR quenching | Identity continuum primitive, separate ER port; ER constructor rejects an NR quenching argument; variable-width ER coordinate checked by repository tests | PASS |
| Observed inputs bypass response/efficiency | Observed continuum, fixed/per-exposure bins, positive/negative observed point locations and endpoint assignment; no response argument on observed-line class | PASS |
| Independent per-component responses | A mixed NR/ER/observed composite has distinct Gaussian widths and efficiencies; its B stays fixed when unrelated signal count provider changes | PASS |
| Monoenergetic areas/probabilities | NR/ER, zero/narrow/broad widths, constant and variable efficiencies, line at/nextafter bin boundary, E=0 threshold, finite-window loss | PASS |
| Lines are not truth quadrature | A policy that raises on truth integration still permits the Gaussian point-line path with unweighted probabilities | PASS |
| Separate line and continuum support | A line outside the continuum range is retained; a continuum callback verifies it is never evaluated outside its own support | PASS |
| Custom ports and matrices | ER-tagged external continuum/line provider and NR kernel line port work; matrix point-lines explicitly rejected because within-bin averages do not determine point response | PASS |
| Compatibility, units, exposure and ownership | Existing observed classes unchanged; fixed exposure rules preserved; domain/type mismatches and invalid line inputs rejected; owned tuples/output copies | PASS |

There are **39 new routing cases**. They include explicit line-area conservation and the expected lost Gaussian probability. Positive sigma assigns half the Gaussian mass on either side of a mean at a bin boundary when the opposite edges are sufficiently far; sigma=0 preserves the documented deterministic `[a,b)` convention with final-edge inclusion. Threshold/tail losses are not restored by renormalization.

`EnergyDomain.NR`, `.ER`, and `.OBSERVED` are machine-readable enum values. They prevent the tested route mismatches and ambiguous string tags. They cannot independently prove that a user's numeric spectrum is physically expressed in its declared unit; external response/spectrum purity and provenance remain caller responsibilities.

### Mixed-background pipeline

The independent composite uses bins `(0,0.2,0.5,1.2)` and X=11:

- NR continuum `2+E` on `[0,3]`, line at E_nr=1.7 with intensity 0.8, q=0.2, sigma=0.04, efficiencies 0.6 and 0.7;
- ER continuum `0.7+E²` on `[0,1.2]`, line at E_er=0.4 with intensity 1.3, identity mean, sigma=0.09, efficiencies 0.8 and 0.5;
- observed continuum `0.3+0.1y` and observed point at 0.6 with intensity 0.9.

Independent truth integrals plus Gaussian line probabilities and observed polynomial/point contributions give:

| Bin | Independent B | Production B |
| --- | ---: | ---: |
| [0,0.2) | 12.141612347452217 | 12.141612347452213 |
| [0.2,0.5) | 36.76629486007172 | 36.76629486007172 |
| [0.5,1.2] | 28.68803361360635 | 28.68803361360635 |

Maximum absolute difference is approximately `3.6e-15` events. The statistical driver receives these same backgrounds when the unrelated signal provider returns either zero or five counts per bin. The signal response is not accidentally applied to B.

## 4. B — SI normalization

### Conversion and frozen coupling convention

From exact SI definitions,

\[
C=\left(\frac{h c_{\rm cm/s}}{2\pi e\,10^9}\right)^2
\simeq3.8937937217185937\times10^{-28},
\]

with h=`6.62607015e-34` J s, e=`1.602176634e-19` J/eV, and c=`29979245800` cm/s. This independently supplies the GeV⁻²-to-cm² conversion and agrees with the existing package constant at relative tolerance `5e-16`. The implementation reuses that constant rather than introducing another conversion.

For explicit proton normalization mass mp and reduced mass `mu_p=mchi*mp/(mchi+mp)`, the frozen free-proton result gives

\[
\sigma_p^{\rm SI}[\mathrm{cm^2}]
=C\frac{\mu_p^2}{\pi}(c_1^p)^2,
\qquad
c_1^p=s_p\frac{\sqrt{\pi\sigma_p/C}}{\mu_p},\quad s_p=\pm1.
\]

For the explicitly selected ratio r=cn/cp,

\[
c_1^0=\frac{(1+r)c_p}{2},\qquad
c_1^1=\frac{(1-r)c_p}{2}.
\]

Every other operator is zero. Equal p/n coupling gives c0=cp,c1=0; proton-only gives c0=c1=cp/2. Negative r preserves the relative sign, and a common sign reversal leaves O1-only rates unchanged. No factor four, weak-scale normalization or full-sum convention was found.

Forty new high-precision signed forward/inverse cases cover zero and tiny/ordinary cross sections, ordinary and highly unequal finite masses. Ten additional cases check O1-only half-isospin coefficients and p/n inversion. Invalid/missing signs, masses, coefficients, ratios and signal configurations are checked separately. Genuine nonrepresentable outputs fail explicitly. Representable subnormal precision and the inherited response/quadrature range are not asserted to be arbitrary precision.

### Coherent absolute cross sections

Synthetic nuclear charges are `h=(A,Z−N)`, with

\[
W_M^{tu}=\frac{2J+1}{4\pi}h_t h_u,
\]

and all other channel entries explicitly zero. Independent contraction gives

\[
P=(Zc_p+Nc_n)^2,
\qquad
\sigma_T=\sigma_p\left(\frac{\mu_T}{\mu_p}\right)^2[Z+rN]^2.
\]

Four synthetic targets/couplings, including A=Z=1 with target mass exactly equal to the supplied mp, pass. The generic differential cross section is integrated over the elastic endpoint at two incoming speeds. The **keV-to-GeV Jacobian `10^-6` and area conversion C** are applied independently to obtain cm². The free-proton total reproduces sigma_p and toy nuclei reproduce the coherent formula. No physical Ge coefficients are used.

### Absolute rates and strength limits

For normalized `f1(beta)=3 beta²/L³`,

\[
\eta(E)=\frac{3}{2L}\max(1-E/E_{\max},0),
\qquad E_{\max}=\frac{2\mu_T^2L^2}{m_T}10^6\ \mathrm{keV}.
\]

Substitution of the proton cross-section relation gives the independent physical rate

\[
R(E)=N_T\frac{\rho_\chi}{m_\chi}c_{\rm cm/s}86400\,
\sigma_p\,10^{-6}\frac{m_T}{2\mu_p^2}[Z+rN]^2\eta(E).
\]

The conversion C cancels exactly against its occurrence in cp²; its correctness is tested separately above. Eight new signed/ratio cases compare complete scalar/array rates at zero, nonzero and inaccessible recoil, at several EFT reference masses and two WIMP spins. The result is independent of numerical `m_N_reference_GeV` for O1; the distinct mp still controls normalization. Shape/empty/owned-output behavior and constructor convention/completeness gates pass.

`cross_section_limit_cm2(mu_upper)` correctly maps strength, not its square root:

\[
\sigma_{p,\rm upper}=\mu_{\rm upper}\sigma_{p,\rm ref}.
\]

The callable coefficients scale as the square root of sigma_ref, while the rate/counts scale linearly. Three full exact-counting pipelines rescale sigma_ref across six orders of magnitude and retain the same physical bound.

| sigma_ref, cm² | Production mu_90 | Production sigma_p,90, cm² |
| ---: | ---: | ---: |
| `2.3e-44` | 478.4530965249441 | `1.1004421220073714e-41` |
| `1.61e-41` | 0.6835044236070633 | `1.1004421220073721e-41` |
| `2.3e-38` | 0.0004784530965249443 | `1.1004421220073718e-41` |

This is a **synthetic normalization result, not a physical sensitivity prediction**. The fixture is A=31,Z=12,mT=23 GeV,J=3/2,mchi=17 GeV,mp=0.91 GeV,r=0.4, rho=0.43, `2.8e25` targets/kg,L=0.006. The detector uses q=0.3, no smearing, efficiencies 0.6 and 0.7, X=37, bins `(0,0.3,0.9,1.5)`, truth `[0,5]`. Background and observations are zero. Independent counts integrate the linear recoil spectrum exactly; `mu_90*sum(S)=−log(0.1)`.

At the smallest reference the independent S is `(0.0009690023366236497,0.0019282697067858996,0.001915289751504033)` and production differs by less than `9e-19` events per bin. The independently calculated physical bound is approximately `1.1004421220073711e-41 cm²`.

The shape-backend test retains `asymptotic_unvalidated` metadata and its own one-sided threshold. Zero exposure gives `no_sensitivity` and maps infinity to infinity. The conversion does not erase LimitResult metadata or upgrade the statistical construction; callers must keep that record.

There are **98 new SI cases**, plus one architecture case, in the final independent test file.

## 5. Architecture, limitations, and findings

The frozen convention document's SHA-256 remains

`a2c57658add45b273ca082d8b01fbab7dd739f8625d6f04a438d619d4cf385cd`.

Git diff confirms no changes to the conversion constant, interactions, kinematics, generic rate engine, statistics, sensitivity or historical notebook. An AST comparison also confirms that every pre-existing function and method body in `detector/response.py` is identical to the accepted base; the change adds domain/point-line methods rather than modifying continuum folding. New background code depends on a structural response protocol and the generic energy enum, not a concrete detector or WIMP model. The SI module is an adapter delegating to the existing rate engine, not a competing rate implementation. Synthetic data remain in tests/examples; no physical Ge or empirical calibration is introduced.

**New BLOCKER/MAJOR/MINOR findings: none.** No new normative-theory contradiction, event-normalization error, routing defect, or ordinary-range numerical failure was identified.

Retained conditions are not silently closed:

1. **Asymptotic shape coverage:** low-count quantitative coverage is not certified; output remains explicitly asymptotic. Exact counting is a distinct construction.
2. **Root coordinate:** shape roots use absolute mu tolerance; tiny strengths require suitable template normalization/tolerances and convergence checks.
3. **Extreme numerical tails/products:** the accepted generic Gaussian tail bound remains. The unchanged remediation stress suite retains its strict failing remote-tail assertion. Neither these additions nor the SI adapter repair upstream erased responses/probabilities, arbitrary extreme products, or ill-conditioned coherent cancellation.

Pure-neutron cp=0 with cn nonzero is outside a finite cn/cp proton-reference parametrization, as explicitly documented; use generic Wilson coefficients for it. Fixed-binned backgrounds cannot be rebinned without more information, and matrix within-bin averages do not determine point-line migration. These are explicit information/model boundaries, not missing physical data or defects.

During new-harness development, two invalid-input assertions mistakenly treated the documented +infinity limit mapping and a finite large coupling ratio as necessarily invalid. They were corrected to honor the declared contract and to probe a genuinely overflowing coefficient configuration instead. Scientific expected values and tolerances were unchanged; the final entire test file was rerun. No production/repository/inherited audit test was changed.

## 6. Execution results and acceptance

All required runs completed. Counts are pytest cases, not a count of individual numerical samples.

| Suite | Exact final result |
| --- | --- |
| Repository | **730 passed in 251.29 s** |
| Unchanged pre-rate, rate, generic audits | **515 passed in 501.40 s** |
| Unchanged detector-remediation extension | **50 passed, 1 failed in 673.42 s**; same documented remote-tail probe |
| New independent routing/SI suite | **138 passed in 128.55 s** |

The boundary failure is `test_local_piecewise_partition[0.009-0.15]`: a weighted tail at about 37.78 sigma has independent mass `5.26582285e-313` but returns zero. It matches the accepted numerical-range observation from CODE_AUDIT_004, with unchanged Gaussian-resolution source. It is retained as a failing assertion, not suppressed, skipped or counted as a pass. The original seven ordinary detector defects remain closed inside the 515-case suite. Test times reflect concurrent audit runs and are not a performance certification.

 The repository has **730 passing cases**, including all 675 inherited cases and 55 new coder cases. The four unchanged original audit files comprise **515 cases**. The unchanged remediation extension has **51 cases**, including its already disclosed strict remote-tail probe. The final new audit has **138 cases: 39 routing, 98 SI, one architecture**.

| Addition / gate | Verdict |
| --- | --- |
| Exact delivered source tree | PASS; commit metadata not independently recovered |
| NR/ER/observed background routing | **PASS** |
| Monoenergetic lines and area/exposure conservation | **PASS** |
| Existing background compatibility and per-component responses | **PASS** |
| Absolute SI conversion / half-isospin coefficients | **PASS** |
| Coherent total cross sections / complete SI rates | **PASS** |
| Physical strength-limit mapping and reference invariance | **PASS** |
| Invalid input, ownership and dependency checks | **PASS** |
| Physical Ge data/calibration/coverage certification | Outside scope; not claimed |
| Overall delivered tree | **PASS WITH CONDITIONS** inherited from the accepted generic foundation |

The verified tree is suitable for subsequent separately sourced physical-input work within those bounds. This certificate does not authenticate the absent commit object, introduce physical inputs, establish experimental calibration, or certify low-count shape coverage. A released commit should preserve the verified tree and its source identity.

## Appendix A. Runnable independent tests

The separately delivered `test_background_si_audit.py` is reproduced verbatim here. Save it in `audit5-work/` next to a candidate checkout named `audit-repo`. No repository-test helper is imported.

```python
"""Independent BACKGROUND_SI_001 audit. Synthetic charges/masses only."""
import ast,json,math
from pathlib import Path
from decimal import Decimal as D,localcontext
from dataclasses import replace,FrozenInstanceError
import mpmath as mp
import numpy as np
import pytest
from scipy.integrate import quad
from wimp_legend import constants
from wimp_legend.energy import EnergyDomain as ED
from wimp_legend.backgrounds import (TruthBackground,MonoenergeticLine,ObservedLineBackground,
    FlatBackground,CallableBackground,BinnedBackground,CompositeBackground,ZeroBackground)
from wimp_legend.detector import (DetectorResponse,ElectronRecoilResponse,ConstantQuenching,
    LindhardQuenching,NoSmearing,GaussianResolution,ConstantEfficiency,CallableEfficiency,
    MatrixResponse,KernelResponse,IntegrationPolicy)
from wimp_legend.si import (si_proton_coefficient_GeV_minus2 as cp_from_sigma,
    si_proton_cross_section_cm2 as sigma_from_cp,si_wilson_coefficients as wilson,SIO1Signal)
from wimp_legend.rates import NREFTRateConfig,differential_cross_section_GeV_minus3 as cross,transition_probability
from wimp_legend.targets import (SourceReference,Isotope,ResponseMetadata,MomentumVariable,
    ResponseKey,PolynomialResponse,NuclearResponseDataset)
from wimp_legend.statistics import PoissonCLsCounting,BinnedPoissonProfile,RootPolicy
from wimp_legend.sensitivity import SensitivityAnalysis
mp.mp.dps=90
CHANNELS=('M','Sigma_double_prime','Sigma_prime','Phi_double_prime','Phi_double_prime_M','Phi_tilde_prime','Delta','Delta_Sigma_prime')
ID='anand2014_prc89_065501_half_isospin_v1'

def conversion():
    p=mp.pi;hc=mp.mpf('6.62607015e-34')*mp.mpf('29979245800')/(2*p*mp.mpf('1.602176634e-10'))
    return hc*hc

def mu(m,M):return mp.mpf(float(m))*mp.mpf(float(M))/(mp.mpf(float(m))+mp.mpf(float(M)))
def coefficient(sigma,m=17.,M=.91,sign=-1):return float(sign*mp.sqrt(mp.pi*mp.mpf(float(sigma))/conversion())/mu(m,M))
def normal_interval(a,b):
    a,b=mp.mpf(float(a)),mp.mpf(float(b));s=mp.sqrt(2)
    if a>=0:return float((mp.erfc(a/s)-mp.erfc(b/s))/2)
    if b<=0:return float((mp.erfc(-b/s)-mp.erfc(-a/s))/2)
    return float((mp.erf(b/s)-mp.erf(a/s))/2)
def weighted(mean,sigma,a,b):
    za,zb=(a-mean)/sigma,(b-mean)/sigma
    p=normal_interval(za,zb);phi=lambda z:math.exp(-z*z/2)/math.sqrt(2*math.pi)
    return (.5+.1*mean)*p+.1*sigma*(phi(za)-phi(zb))

def nucleus(A=31,Z=12,M=23.,J=1.5):
    source=SourceReference('Independent synthetic coherent charges','audit-only algebra','synthetic')
    iso=Isotope('X',A,Z,M,J,dict.fromkeys(('identity','mass_GeV','spin'),source))
    metadata=ResponseMetadata('synthetic coherent','0 1'.split(),'half-isospin','raw W','dimensionless',MomentumVariable('q2_GeV2'),source,ID)
    h=(A,2*Z-A);responses={}
    for ch in CHANNELS:
        for t in (0,1):
            for u in (0,1):
                v=(2*J+1)/(4*math.pi)*h[t]*h[u] if ch=='M' else 0.
                responses[ResponseKey(ch,t,u)]=PolynomialResponse((v,),0,None,source)
    return NuclearResponseDataset('independent-synthetic-SI',iso,metadata,responses)

def signal(sigma=2.3e-42,ratio=.4,sign=-1,ds=None,mref=1.7,**changes):
    kw=dict(sigma_p_ref_cm2=sigma,m_chi_GeV=17.,m_proton_GeV=.91,fn_over_fp=ratio,proton_sign=sign,
        dataset=nucleus() if ds is None else ds,rho_chi_GeV_cm3=.43,targets_per_kg=2.8e25,
        j_chi=1.,config=NREFTRateConfig(mref),speed_pdf=lambda b:3*(b/.006)**2/.006,beta_max=.006)
    kw.update(changes);return SIO1Signal(**kw)

def rate_ref(E,sigma=2.3e-42,ratio=.4,A=31,Z=12,M=23.,m=17.,p=.91):
    mt=float(mu(m,M));mp_=float(mu(m,p));L=.006
    endpoint=2*mt*mt*L*L/M*1e6
    K=2.8e25*.43/m*29979245800*86400*sigma*1e-6*M/(2*mp_**2)*(Z+ratio*(A-Z))**2*3/(2*L)
    return K*np.maximum(1-np.asarray(E)/endpoint,0)

# A. Explicit routing / point measures.
@pytest.mark.parametrize('domain,q',[(ED.NR,.2),(ED.NR,.65),(ED.ER,1.)])
@pytest.mark.parametrize('X',[0.,1.,11.])
def test_continuum_route_jacobian(domain,q,X):
    response=DetectorResponse(ConstantQuenching(q),NoSmearing()) if domain is ED.NR else ElectronRecoilResponse(NoSmearing())
    b=TruthBackground(domain,response,spectrum=lambda E:2+E,spectrum_energy_range_keV=(0,3))
    bins=np.array([0.,.2,.5,1.2]);l=np.clip(bins/q,0,3);F=lambda e:2*e+e*e/2
    expected=X*np.diff(F(l))
    np.testing.assert_allclose(b.expected_counts(bins,X),expected,rtol=3e-14,atol=0)

@pytest.mark.parametrize('domain,q',[(ED.NR,.23),(ED.ER,1.)])
@pytest.mark.parametrize('sigma',[0.,.03,.4])
def test_truth_line_exact_area_and_efficiency(domain,q,sigma):
    E=.8;I=1.7;X=9.;bins=np.array([-.2,.1,.4,.8,1.4]);mean=q*E
    # Truth efficiency depends on truth energy, observed efficiency on y.
    kw=dict(truth_efficiency=CallableEfficiency(lambda e:.2+.1*e),observed_efficiency=CallableEfficiency(lambda y:.5+.1*y))
    response=DetectorResponse(ConstantQuenching(q),GaussianResolution(sigma),**kw) if domain is ED.NR else ElectronRecoilResponse(GaussianResolution(sigma),**kw)
    b=TruthBackground(domain,response,lines=[MonoenergeticLine(E,I,domain)])
    ref=np.zeros(len(bins)-1)
    if sigma:
        ref=np.array([weighted(mean,sigma,a,c) for a,c in zip(bins[:-1],bins[1:])])
    elif bins[0]<=mean<=bins[-1]:
        j=min(np.searchsorted(bins,mean,side='right')-1,len(ref)-1);ref[j]=.5+.1*mean
    ref*=X*I*(.2+.1*E)
    np.testing.assert_allclose(b.expected_counts(bins,X),ref,rtol=3e-8,atol=1e-13)

@pytest.mark.parametrize('delta',[-1,0,1])
def test_line_threshold_nextafter(delta):
    edge=.4;E=np.nextafter(edge,0.) if delta<0 else np.nextafter(edge,1.) if delta>0 else edge
    bins=[0,.4,.8]
    for width in [0.,1e-7]:
        b=TruthBackground(ED.ER,ElectronRecoilResponse(GaussianResolution(width)),lines=[MonoenergeticLine(E,2.,ED.ER)])
        ref=np.array([normal_interval((a-E)/width,(c-E)/width) for a,c in zip(bins[:-1],bins[1:])]) if width else np.array([1,0] if delta<0 else [0,1])
        np.testing.assert_allclose(b.expected_counts(bins,3.),6*ref,rtol=2e-13,atol=0)

@pytest.mark.parametrize('width',[.03,.4])
def test_line_conservation_partial_window_and_no_delta_quadrature(width):
    # Disallow any TRUTH quadrature; only Gaussian bin probabilities are needed.
    class NoTruthIntegration(IntegrationPolicy):
        def integrate(self,*a,**k):raise AssertionError('line must not use truth quadrature')
    d=DetectorResponse(ConstantQuenching(.2),GaussianResolution(width),policy=NoTruthIntegration())
    b=TruthBackground(ED.NR,d,lines=[MonoenergeticLine(2.,3.,ED.NR)])
    bins=[-3,.4,4]
    ref=15*np.array([normal_interval((-3-.4)/width,0),normal_interval(0,(4-.4)/width)])
    np.testing.assert_allclose(b.expected_counts(bins,5),ref,rtol=3e-14)
    assert sum(b.expected_counts([.4,.41],5))<7.5
    assert sum(b.expected_counts(bins,5))==pytest.approx(sum(ref),rel=3e-14)

def test_line_outside_continuum_support_and_kernel_port():
    def spec(E):
        assert 0<=E<=1;return 1.
    d=DetectorResponse(ConstantQuenching(.2),NoSmearing())
    b=TruthBackground(ED.NR,d,spectrum=spec,spectrum_energy_range_keV=(0,1),lines=[MonoenergeticLine(4.,2.,ED.NR)])
    np.testing.assert_allclose(b.expected_counts([0,.4,1.],3),[3.,6.],rtol=2e-14)
    k=KernelResponse(lambda E,bins:[.1,.7],True)
    c=TruthBackground(ED.NR,k,lines=[MonoenergeticLine(2.,3.,ED.NR)])
    np.testing.assert_allclose(c.expected_counts([0,1,2],4),[1.2,8.4],rtol=2e-14)
    with pytest.raises(TypeError):TruthBackground(ED.NR,MatrixResponse([0,1],[0,1],[[1]],True,'synthetic'),lines=[MonoenergeticLine(.5,1,ED.NR)])

def test_observed_bypass_and_legacy_exposure_contracts():
    bins=[-.2,0,.5,1.];X=7.
    lines=[MonoenergeticLine(x,1.,ED.OBSERVED) for x in [-.2,0,.5,1,1.1]]
    b=ObservedLineBackground(lines)
    np.testing.assert_array_equal(b.expected_counts(bins,X),[7,7,14])
    components=CompositeBackground([FlatBackground(.3),CallableBackground(lambda y:.2+.1*y),b])
    F=lambda y:.5*y+.05*y*y
    ref=X*np.diff(F(np.array(bins)))+np.array([7,7,14])
    np.testing.assert_allclose(components.expected_counts(bins,X),ref,rtol=3e-15)
    rates=BinnedBackground(bins,[1,2,3],'per_unit_exposure')
    fixed=BinnedBackground(bins,[1,2,3],'fixed_expected_counts',7)
    np.testing.assert_array_equal(rates.expected_counts(bins,14),[14,28,42])
    np.testing.assert_array_equal(fixed.expected_counts(bins,7),[1,2,3])
    with pytest.raises(ValueError):fixed.expected_counts(bins,14)


def test_composite_independent_responses_and_signal_response_isolation():
    bins=np.array([0.,.2,.5,1.2]);X=11.
    nr=DetectorResponse(ConstantQuenching(.2),GaussianResolution(.04),truth_efficiency=ConstantEfficiency(.6),observed_efficiency=ConstantEfficiency(.7))
    er=ElectronRecoilResponse(GaussianResolution(.09),truth_efficiency=ConstantEfficiency(.8),observed_efficiency=ConstantEfficiency(.5))
    components=[TruthBackground(ED.NR,nr,spectrum=lambda E:2+E,spectrum_energy_range_keV=(0,3),lines=[MonoenergeticLine(1.7,.8,ED.NR)]),
        TruthBackground(ED.ER,er,spectrum=lambda E:.7+E*E,spectrum_energy_range_keV=(0,1.2),lines=[MonoenergeticLine(.4,1.3,ED.ER)]),
        CallableBackground(lambda y:.3+.1*y),ObservedLineBackground([MonoenergeticLine(.6,.9,ED.OBSERVED)])]
    expected=[]
    for a,b in zip(bins[:-1],bins[1:]):
        pnr=lambda E:normal_interval((a-.2*E)/.04,(b-.2*E)/.04)
        per=lambda E:normal_interval((a-E)/.09,(b-E)/.09)
        n=.42*(quad(lambda E:(2+E)*pnr(E),0,3,epsabs=1e-12)[0]+.8*pnr(1.7))
        e=.4*(quad(lambda E:(.7+E*E)*per(E),0,1.2,epsabs=1e-12)[0]+1.3*per(.4))
        o=.3*(b-a)+.05*(b*b-a*a)+(.9 if a<=.6<b else 0)
        expected.append(X*(n+e+o))
    background=CompositeBackground(components)
    # Deliberately unrelated signal response. Background must retain its OWN responses.
    class ArbitrarySignal:
        def __init__(self,v):self.v=v
        def expected_counts(self,*a,**k):return np.full(3,self.v)
    for value in [0.,5.]:
        a=SensitivityAnalysis(ArbitrarySignal(value),background,PoissonCLsCounting())
        r=a.observed_limit(lambda E:1,bins,X,(0,3),[1,2,3])
        np.testing.assert_allclose(r.background_counts,expected,rtol=2e-10,atol=1e-12)
        np.testing.assert_array_equal(r.signal_counts,[value]*3)
    Path(__file__).with_name('background-pipeline.json').write_text(json.dumps(dict(expected=expected,actual=list(r.background_counts)),indent=2))

@pytest.mark.parametrize('bad',[True,'NR','keV_nr',0,np.array('NR'),ED.OBSERVED])
def test_domain_mismatch_rejected(bad):
    with pytest.raises((ValueError,TypeError)):TruthBackground(bad,DetectorResponse(ConstantQuenching(.2),NoSmearing()),lines=[MonoenergeticLine(1,1,ED.NR)])
    with pytest.raises((ValueError,TypeError)):MonoenergeticLine(1,1,bad) if bad is not ED.OBSERVED else TruthBackground(ED.NR,ElectronRecoilResponse(NoSmearing()),lines=[MonoenergeticLine(1,1,ED.NR)])

@pytest.mark.parametrize('bad',[-1.,np.nan,np.inf,True,'1',1j,[1]])
def test_bad_lines(bad):
    with pytest.raises((ValueError,TypeError)):MonoenergeticLine(bad,1,ED.NR)
    with pytest.raises((ValueError,TypeError)):MonoenergeticLine(1,bad,ED.ER)
    with pytest.raises((ValueError,TypeError)):ObservedLineBackground([MonoenergeticLine(1,1,ED.NR)])


def test_custom_ER_port_and_owned_line_inputs():
    class CustomER:
        input_energy_domain=ED.ER
        def expected_counts(self,spectrum,bins,X,interval,**kwargs):return np.array([.3,.4])*X
        def expected_line_counts(self,E,I,bins,X):return X*I*np.array([.2,.6])
    lines=[MonoenergeticLine(.5,2,ED.ER)]
    b=TruthBackground(ED.ER,CustomER(),spectrum=lambda E:1,spectrum_energy_range_keV=[0,1],lines=lines)
    lines.clear();a=b.expected_counts([0,.4,1],3);np.testing.assert_allclose(a,[2.1,4.8],rtol=2e-15)
    a[:]=9;np.testing.assert_allclose(b.expected_counts([0,.4,1],3),[2.1,4.8],rtol=2e-15)
    with pytest.raises(FrozenInstanceError):b.energy_domain=ED.NR
    with pytest.raises(TypeError):ElectronRecoilResponse(NoSmearing(),quenching=ConstantQuenching(.2))

# B. SI normalization.
def test_SI_defining_constants_and_no_hidden_factor():
    assert constants.GEV_MINUS2_TO_CM2==pytest.approx(float(conversion()),rel=5e-16)
    assert float(conversion())==pytest.approx(3.8937937217185937e-28,rel=5e-16)

@pytest.mark.parametrize('m,p',[(.02,.91),(17.,.91),(1e7,1.2),(1e-100,1e100),(1e100,1e-100)])
@pytest.mark.parametrize('sigma',[0.,1e-300,2.3e-42,1e-20])
@pytest.mark.parametrize('sign',[-1,1])
def test_forward_inverse_scalar_math(m,p,sigma,sign):
    ref=coefficient(sigma,m,p,sign)
    got=cp_from_sigma(sigma,m_chi_GeV=m,m_proton_GeV=p,proton_sign=sign)
    assert got==pytest.approx(ref,rel=1e-14,abs=0)
    back=sigma_from_cp(got,m_chi_GeV=m,m_proton_GeV=p)
    assert back==pytest.approx(sigma,rel=2e-14,abs=0)
    assert sigma_from_cp(-got,m_chi_GeV=m,m_proton_GeV=p)==back

@pytest.mark.parametrize('r',[1.,0.,-1.,.4,-.7])
@pytest.mark.parametrize('sign',[-1,1])
def test_half_isospin_O1_only(r,sign):
    sigma=2.3e-42;cp=coefficient(sigma,sign=sign)
    w=wilson(sigma,m_chi_GeV=17,m_proton_GeV=.91,fn_over_fp=r,proton_sign=sign)
    assert w.convention_id==ID
    np.testing.assert_allclose(w.values[1],[(1+r)*cp/2,(1-r)*cp/2],rtol=3e-15,atol=0)
    assert all(pair==(0.,0.) for op,pair in w.values.items() if op!=1)
    pn=w.to_proton_neutron()[1];np.testing.assert_allclose(pn,[cp,r*cp],rtol=3e-15,atol=abs(cp)*1e-15)

@pytest.mark.parametrize('A,Z,M,J,r',[(1,1,.91,.5,0.),(31,12,23.,1.5,.4),(22,9,20.,0.,-.3),(20,10,18.,.5,1.)])
def test_coherent_absolute_total_cross_sections(A,Z,M,J,r):
    sig=2.3e-42;ds=nucleus(A,Z,M,J);w=wilson(sig,m_chi_GeV=17.,m_proton_GeV=.91,fn_over_fp=r,proton_sign=-1)
    cp=coefficient(sig);amplitude=cp*cp*(Z+r*(A-Z))**2
    got=transition_probability(ds,w,q2_GeV2=0,V=0,j_chi=1,config=NREFTRateConfig(1.7))
    assert got==pytest.approx(amplitude,rel=3e-14,abs=0)
    expected=sig*float(mu(17,M)/mu(17,.91))**2*(Z+r*(A-Z))**2
    for beta in [.002,.009]:
        end=2*float(mu(17,M))**2*beta**2/M*1e6
        natural=quad(lambda E:cross(E,beta,m_chi_GeV=17.,dataset=ds,coefficients=w,j_chi=1.,config=NREFTRateConfig(1.7)),0,end,epsabs=0,epsrel=1e-11)[0]
        actual=natural*1e-6*float(conversion())
        assert actual==pytest.approx(expected,rel=4e-14,abs=0)

@pytest.mark.parametrize('r',[0.,.4,1.,-.3])
@pytest.mark.parametrize('sign',[-1,1])
def test_absolute_SI_recoil_rates(r,sign):
    s=signal(ratio=r,sign=sign);E=np.array([[0.,.1,1.],[5.,20.,1e8]])
    np.testing.assert_allclose(s(E),rate_ref(E,ratio=r),rtol=4e-13,atol=0)
    assert s(np.empty((0,2))).shape==(0,2)
    assert isinstance(s(np.array(1.)),float)
    for mref in [.3,3.]:np.testing.assert_allclose(signal(ratio=r,sign=sign,mref=mref)(E),rate_ref(E,ratio=r),rtol=4e-13,atol=0)
    s2=signal(ratio=r,sign=-sign,j_chi=0.)
    np.testing.assert_allclose(s2(E),s(E),rtol=4e-14,atol=0)

@pytest.mark.parametrize('factor',[.01,7.,1e4])
def test_sigma_scaling_and_exact_count_limit(factor):
    sig=2.3e-42*factor;s=signal(sigma=sig)
    bins=np.array([0.,.3,.9,1.5]);q=.3;E=bins/q
    K=float(rate_ref(0,sigma=sig));end=2*float(mu(17,23))**2*.006**2/23*1e6
    F=lambda x:K*(x-x*x/(2*end))
    expected=37*.42*np.diff(F(E))
    analysis=SensitivityAnalysis(DetectorResponse(ConstantQuenching(q),NoSmearing(),truth_efficiency=ConstantEfficiency(.6),observed_efficiency=ConstantEfficiency(.7)),ZeroBackground(),PoissonCLsCounting())
    result=analysis.observed_limit(s,bins,37.,(0,5),[0,0,0])
    np.testing.assert_allclose(result.signal_counts,expected,rtol=8e-13)
    assert result.limit.mu_upper*sum(expected)==pytest.approx(-math.log(.1),rel=8e-13)
    actual=s.cross_section_limit_cm2(result.limit.mu_upper)
    ref=-math.log(.1)/sum(expected)*sig
    assert actual==pytest.approx(ref,rel=8e-13,abs=0)
    cp=coefficient(sig);assert s.coefficients.to_proton_neutron()[1][0]==pytest.approx(cp,rel=2e-15)
    assert signal().cross_section_limit_cm2(0)==0
    Path(__file__).with_name('si-pipeline-'+str(factor)+'.json').write_text(json.dumps(dict(sigma_ref=sig,S_reference=expected.tolist(),S_actual=result.signal_counts,mu_upper=result.limit.mu_upper,sigma_limit_reference=ref,sigma_limit_actual=actual,method=result.limit.method),indent=2))


def test_shape_metadata_zero_sensitivity_and_mutation():
    s=signal();bins=[0,.5,1.]
    d=DetectorResponse(ConstantQuenching(.3),NoSmearing())
    a=SensitivityAnalysis(d,ZeroBackground(),BinnedPoissonProfile(RootPolicy(xtol=1e-20,rtol=1e-12)))
    result=a.observed_limit(s,bins,10.,(0,4),[0,0])
    assert result.limit.status=='asymptotic_unvalidated'
    z=float(mp.sqrt(2)*mp.erfinv(mp.mpf('.8')))
    assert result.limit.mu_upper*sum(result.signal_counts)==pytest.approx(z*z/2,rel=3e-12)
    empty=SensitivityAnalysis(d,ZeroBackground(),PoissonCLsCounting()).expected_limit(s,bins,0.,(0,4))
    assert empty.limit.status=='no_sensitivity' and s.cross_section_limit_cm2(empty.limit.mu_upper)==np.inf
    E=np.arange(6.).reshape(2,3)[:,::2];original=E.copy();out=s(E);saved=out.copy();E[:]=99
    np.testing.assert_array_equal(out,saved)
    with pytest.raises(FrozenInstanceError):s.sigma_p_ref_cm2=1.
    assert signal().convention_id==ID

@pytest.mark.parametrize('bad',[-1.,0.,np.nan,np.inf,True,'1',1j,[1]])
def test_invalid_masses(bad):
    for field in ['m_chi_GeV','m_proton_GeV']:
        kw=dict(m_chi_GeV=17.,m_proton_GeV=.91);kw[field]=bad
        with pytest.raises((ValueError,TypeError)):cp_from_sigma(1e-42,proton_sign=1,**kw)
        with pytest.raises((ValueError,TypeError)):sigma_from_cp(1e-7,**kw)
        with pytest.raises((ValueError,TypeError)):signal(**{field:bad})

@pytest.mark.parametrize('bad',[-1.,np.nan,np.inf,True,'1',1j,[1]])
def test_invalid_sigma_or_coefficient(bad):
    with pytest.raises((ValueError,TypeError)):cp_from_sigma(bad,m_chi_GeV=17,m_proton_GeV=.91,proton_sign=1)
    if bad!=-1.:
        with pytest.raises((ValueError,TypeError)):sigma_from_cp(bad,m_chi_GeV=17,m_proton_GeV=.91)
    with pytest.raises((ValueError,TypeError)):signal(sigma=bad)
    if bad==np.inf:
        assert signal().cross_section_limit_cm2(bad)==np.inf
    else:
        with pytest.raises((ValueError,TypeError)):signal().cross_section_limit_cm2(bad)

@pytest.mark.parametrize('bad',[0,2,-2,True,1.,-1.,'1',np.nan])
def test_bad_or_missing_sign(bad):
    with pytest.raises((ValueError,TypeError)):cp_from_sigma(1e-42,m_chi_GeV=17,m_proton_GeV=.91,proton_sign=bad)
    with pytest.raises(TypeError):cp_from_sigma(1e-42,m_chi_GeV=17,m_proton_GeV=.91)
    with pytest.raises(TypeError):wilson(1e-42,m_chi_GeV=17,m_proton_GeV=.91,proton_sign=1)

@pytest.mark.parametrize('bad',[np.nan,np.inf,True,1j,'1',[1]])
def test_invalid_ratio_and_zero_reference(bad):
    with pytest.raises((ValueError,TypeError)):signal(ratio=bad)
    with pytest.raises(ValueError):signal(sigma=0.)
    assert wilson(0.,m_chi_GeV=17,m_proton_GeV=.91,fn_over_fp=.4,proton_sign=-1).values[1]==(0.,0.)


def test_SI_range_failures_are_explicit():
    for cp,m,p in [(1e308,1e308,1e308),(1e-300,1e-300,1e-300)]:
        with pytest.raises((ValueError,OverflowError,FloatingPointError)):sigma_from_cp(cp,m_chi_GeV=m,m_proton_GeV=p)
    with pytest.raises(ValueError):cp_from_sigma(1e308,m_chi_GeV=1e-308,m_proton_GeV=1e-308,proton_sign=1)
    with pytest.raises(ValueError):wilson(1e280,m_chi_GeV=17,m_proton_GeV=.91,fn_over_fp=1e308,proton_sign=1)


def test_architecture_protected_sources_and_purity_contract():
    root=Path(__file__).resolve().parents[1]/'audit-repo/src/wimp_legend'
    for name in ['backgrounds/routing.py','energy.py']:
        for node in ast.walk(ast.parse((root/name).read_text())):
            mods=[n.name for n in node.names] if isinstance(node,ast.Import) else [node.module or ''] if isinstance(node,ast.ImportFrom) else []
            for mod in mods:assert not set(mod.split('.'))&{'detector','si','rates','targets','interactions','halo','statistics','sensitivity'}

@pytest.mark.parametrize('sigma',[0.,.1])
def test_NR_line_exact_threshold_zero_and_boundary(sigma):
    response=DetectorResponse(ConstantQuenching(.2),GaussianResolution(sigma))
    for E in [0.,2.]:
        mean=.2*E;bins=[0,.4,.8]
        actual=response.expected_line_counts(E,3.,bins,5.)
        if sigma:ref=15*np.array([normal_interval((a-mean)/sigma,(b-mean)/sigma) for a,b in zip(bins[:-1],bins[1:])])
        else:ref=np.array([15.,0.]) if E==0 else np.array([0.,15.])
        np.testing.assert_allclose(actual,ref,rtol=3e-14,atol=0)
        np.testing.assert_array_equal(response.expected_line_counts(E,0.,bins,5.),[0,0])
        np.testing.assert_array_equal(response.expected_line_counts(E,3.,bins,0.),[0,0])


def test_SI_constructor_gate_and_input_ownership():
    s=signal()
    mismatched=replace(s.dataset,metadata=replace(s.dataset.metadata,convention_id=ID+' '))
    with pytest.raises(ValueError):signal(ds=mismatched)
    responses=dict(s.dataset.responses);responses.pop(ResponseKey('Delta',1,0))
    with pytest.raises(KeyError):signal(ds=replace(s.dataset,responses=responses))
    points=np.array([.003]);s=signal(breakpoints=points);points[:]=.9
    assert s.breakpoints==(.003,)
    assert all(op==1 or pair==(0.,0.) for op,pair in s.coefficients.values.items())
```

## Appendix B. Reference note frozen before implementation/regression inspection

```text
# BACKGROUND_SI_001 references frozen before new implementation / regression inspection
Source-lock targets from handoff: base b48fc6d040ec30dc3c008a981773bba16bff300f, commit9d2468e875792cd19210904e401d383d242d6cc5, tree2da8deb5cac7c9f07c73722f0288e1ced38a5c57. Verify actual base and patch tree. No code edits beyond applying delivered source.

Backgrounds: continuum r in events/kg/day/keV_truth -> X integral_truth r(E)*eps_truth(E)*p_j(E)dE. NR mean=q(E)E; ER identity mean=E (explicit external ER model may differ). Observed spectrum/count already accepted, never folded again. NR constant q maps density r(y/q)/q; ER has identity Jacobian. For polynomial r(E)=a+bE primitive aE+bE²/2; deterministic NR bin integrals primitive(b/q)-primitive(a/q) intersect explicit truth range, ER q1. Each component owns its response independent of signal response.
Line intensity I is integrated events/kg/day: counts_j=X*I*eps_truth(Eline)*integral_bin Normal(y;mean,sigma)*eps_obs(y)dy. NO inverse quenching Jacobian in line area. Normal probability from independent erfc interval. At sigma>0 mean on bin boundary gives half probability on each side (finite-window tails lost); sigma0 deterministic documented [a,b) finaledge inclusion. Lines can be outside continuum range; line area not integrated by delta quadrature. Observed points direct binning exposure*I, and observed fixedcounts preserve reference exposure. Positive sigma small/broad, threshold nextafter/zero-width, zero weights, independent components, truth/observed efficiencies once. Matrix within-bin average insufficient for unknown point response; reject lines without explicit line port.

SI: exact h=6.62607015e-34 J*s,e=1.602176634e-19 J/eV,c=29979245800 cm/s. C=[h*c/(2*pi*e*1e9)]² ~3.8937937217185937e-28 cm² per GeV^-2. mu_p=mchi*mp/(mchi+mp). sigma_p_cm²=C*mu_p²*cp²/pi. Inverse cp=sign*sqrt(pi*sigma_p/C)/mu_p. cp has GeV^-2, sign must explicit +/-1. cn=r*cp; canonical half isospin c0=cp(1+r)/2,c1=cp(1-r)/2; all other operators zero. cp=cn -> c0=cp,c1=0. neutronzero -> c0=c1=cp/2. Negative r and global sign preserve relative physics. Pure neutron nonzero with cp0 impossible using finite r; generic API appropriate. sigma0 -> cp0 but explicit requiredsign validation still expected. Forward/inverse extremes use high precision and explicit numerical domain checks, no hidden mass, weak-scale, or factor4.

Synthetic coherent response: h=(A,Z-N), W_Mtu=(2J+1)h_t h_u/(4*pi), other28 entries explicit0. P=cp²[Z+rN]². ds/dEGeV=mT P/(2*pi beta²); integrate EkeV0..2muT² beta²/mT*1e6, multiply1e-6*C gives sigmaT=sigmap*(muT/mup)²[Z+rN]². A1Z1 toy targetmass=explicitmp gives exactlysigmap. No nuclear Ge data.

SI fullrate for normalized f=3beta²/L³, muT=mchi*mT/(mchi+mT), Emax=2muT²L²/mT*1e6: eta(E)=3/(2L)*max(1-E/Emax,0). R(E)=NT*rho/mchi*c*86400*sigmap*1e-6*mT/(2mup²)*[Z+rN]²*eta(E), events/kg/day/keV. This uses sigma in cm² and cancels C once against cp. Coefficients sqrt(sigma) and signal proportional sigma; O1 independent EFTreference mN and WIMPspin. Proton normalization mass not EFT mN.

Configuration for independent SI/pipeline: toy A31Z12N19,mT23,J1.5,mchi17,mp.91,r=.4,sign-1,rho.43,NT2.8e25,L.006,sigmaref2.3e-42,mref1.7,jchi1. Signal normalization K from exact SI. Detector q.3 identitysmear, efficiencytruth.6,obs.7,bins[0,.3,.9,1.5],truth0..5,X37; expected S bins X*.42*K*sigmaref*integral(1-E/Emax)dE. Background0 integerobs0: exact CLs mu90*sumS=-ln(.1), sigma90=mu90*sigmaref; rescaling ref by factors.01,7,1e4 preserves physical upper bound. Shapezeroobs z90²/(2sumS) explicitlyapprox; preserve status. Infinite no-sensitivity maps infinity not finite certificate.

Independent background pipeline: qNR.2, rNR(E)=2+E on0..3; ER rate .7+E² on0..1.2 identitymean; separate widths/eff; lineNR E1.7 I.8 ER E.4 I1.3 and observedline.6 I.9. X11 bins[0,.2,.5,1.2]; observedcontinuum .3+.1y alreadyobserved. S/B counts independent sum, change unrelated signalresponse must not alterbackground. Use highprecisionerfc/Gaussian integrals and independent observable-density reference where helpful.
```

## Appendix C. Evidence identity and reproduction

| Evidence | SHA-256 |
| --- | --- |
| `upload/BACKGROUND_SI_001.diff` | `f2f237d0d94cbd2723d7d283bfa8982ca39733160d7247342d1deb06db6c72c3` |
| `upload/BACKGROUND_SI_001_TESTER_HANDOFF.md` | `a69c47f7d4d962c23117aec722ad1c07b75b5c722ca5c0977f860ac45096754f` |
| `upload/BACKGROUND_SI_001_IMPLEMENTATION_REPORT.md` | `29a28f60c1c3e544d386059e7f4b96471a096a3a1226a9423330516e57af4c7e` |
| `audit5-work/test_background_si_audit.py` | `334c363e20464c41892037c8705057ad6cb445c731d9b2e6bc5b6f49602a2ac8` |
| `audit5-work/independent_references.md` | `1c9d7d3486245aaaa6f9a70871bf1f3f0d96fecce7a2ef16f609211b5563636d` |
| `audit4-work/test_independent.py` | `1f1827d261b49a0928db7eec5c33c0bc5d518d976d107f245fd4de6639ac23e6` |
| `audit4-work/test_falsifications.py` | `99438a0db818631638f99d80b7282c8f6ed8d56505443abf4057b4e0005de6cb` |
| `audit4-work/test_rate_audit.py` | `a90623aca9133d9293181cad00c4db5dbbb5fbc0f72cc72793612c172ba24aae` |
| `audit4-work/test_generic_audit.py` | `1705ebe6e24403a54573d3dd6b0ba8724c681bbc29a1e684a93e8bdc1fc05aec` |
| `audit4-work/test_remediation.py` | `99110bc9daf4165844d3439e3f46a6ce85e7acfb6b9086b912a13f4813dd1847` |
| `audit-repo/docs/NREFT_CONVENTIONS.md` | `a2c57658add45b273ca082d8b01fbab7dd739f8625d6f04a438d619d4cf385cd` |
| `audit-repo/src/wimp_legend/si.py` | `e4155a29c2174955821fb87d28f4f50bd37a96f029669e7f0efc25208150aca9` |
| `audit-repo/src/wimp_legend/backgrounds/routing.py` | `5fd754a854ff7edbc5561f5b7b20ec09916f6cf0cb326eebf3845ac90c02fca9` |
| `audit-repo/src/wimp_legend/detector/electron.py` | `5074062e4caae5625887dd52121858ed23e8623f411e5c16a1476b7eb47895ca` |
| `audit-repo/src/wimp_legend/detector/response.py` | `543dda6e1608c612e7b96fa51d768b814d931b8cd9c9308aeec50a0bea677ac0` |

Source reconstruction with the attached diff in `upload/`:

```bash
git clone --branch refactor --single-branch https://github.com/tiagocamison/WIMP-LEGEND.git audit-repo
git -C audit-repo switch --detach b48fc6d040ec30dc3c008a981773bba16bff300f
git -C audit-repo apply --check ../upload/BACKGROUND_SI_001.diff
git -C audit-repo apply --index ../upload/BACKGROUND_SI_001.diff
git -C audit-repo write-tree
# Expected tree: 2da8deb5cac7c9f07c73722f0288e1ced38a5c57
python -m venv --system-site-packages audit-venv
audit-venv/bin/python -m pip install -e './audit-repo[test]' mpmath==1.4.1
audit-venv/bin/python -m pytest -q -W error audit5-work/test_background_si_audit.py
```

Run `../audit-venv/bin/python -m pytest -q -W error` from the checkout for the 730-case repository suite. For closest reproduction pin Python/NumPy/SciPy/pytest to Section 2's versions.

The four old audit files are available verbatim from CODE_AUDIT_001/002/003 Appendices A; the fifth is CODE_AUDIT_004 Appendix A. Keep their recorded hashes and assertions unchanged. Their filesystem-based architecture checks expect `audit-repo` beside their test directory. The boundary file's one strict subnormal-tail failure is expected under the retained documented condition, not an instruction to change it.

Final new-suite output:

```text
........................................................................ [ 52%]
..................................................................       [100%]
138 passed in 128.55s (0:02:08)
```

Unchanged boundary-suite failure summary:

```text
=========================== short test summary info ============================
FAILED audit4-work/test_remediation.py::test_local_piecewise_partition[0.009-0.15]
1 failed, 50 passed in 673.42s (0:11:13)
```

background-pipeline.json:

```json
{
  "expected": [
    12.141612347452217,
    36.76629486007172,
    28.68803361360635
  ],
  "actual": [
    12.141612347452213,
    36.76629486007172,
    28.68803361360635
  ]
}
```

si-pipeline-0.01.json:

```json
{
  "sigma_ref": 2.3e-44,
  "S_reference": [
    0.0009690023366236497,
    0.0019282697067858996,
    0.001915289751504033
  ],
  "S_actual": [
    0.0009690023366236494,
    0.0019282697067858992,
    0.0019152897515040323
  ],
  "mu_upper": 478.4530965249441,
  "sigma_limit_reference": 1.1004421220073711e-41,
  "sigma_limit_actual": 1.1004421220073714e-41,
  "method": "poisson_cls_counting"
}
```

si-pipeline-7.0.json:

```json
{
  "sigma_ref": 1.61e-41,
  "S_reference": [
    0.6783016356365549,
    1.3497887947501297,
    1.3407028260528235
  ],
  "S_actual": [
    0.6783016356365544,
    1.3497887947501288,
    1.3407028260528222
  ],
  "mu_upper": 0.6835044236070633,
  "sigma_limit_reference": 1.100442122007371e-41,
  "sigma_limit_actual": 1.1004421220073721e-41,
  "method": "poisson_cls_counting"
}
```

si-pipeline-10000.0.json:

```json
{
  "sigma_ref": 2.3e-38,
  "S_reference": [
    969.00233662365,
    1928.2697067858999,
    1915.289751504034
  ],
  "S_actual": [
    969.0023366236495,
    1928.2697067858985,
    1915.2897515040318
  ],
  "mu_upper": 0.0004784530965249443,
  "sigma_limit_reference": 1.1004421220073706e-41,
  "sigma_limit_actual": 1.1004421220073718e-41,
  "method": "poisson_cls_counting"
}
```
