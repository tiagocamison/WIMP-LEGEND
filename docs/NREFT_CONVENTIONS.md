# WIMP-LEGEND NREFT conventions

**Status:** FROZEN — theoretical convention v1  
**Date frozen:** 2026-09-13  
**Project convention ID:** `anand2014_prc89_065501_half_isospin_v1`

This document is the canonical theory/convention reference for the WIMP-LEGEND nonrelativistic effective-field-theory (NREFT) implementation. Code, tests, future Astra runs, and nuclear-response datasets should follow this document unless the convention is deliberately revised and the convention ID is changed.

The purpose is to prevent silent convention mismatches between Wilson coefficients, WIMP response functions, nuclear responses, cross sections, and rates.

---

## 1. Canonical sources and authority

### Primary theory source

N. Anand, A. L. Fitzpatrick, and W. C. Haxton,  
**“Weakly interacting massive particle-nucleus elastic scattering response,”**  
*Physical Review C* **89**, 065501 (2014),  
DOI: `10.1103/PhysRevC.89.065501`.

This is the authoritative source for the NREFT operator basis, isospin definitions, WIMP response functions \(R_k^{\tau\tau'}\), nuclear responses \(W_k^{\tau\tau'}\), transition-probability normalization, differential cross section, and rate normalization used here. Relevant equations are especially Eqs. (10)–(18), (37)–(41), (50)–(53), and Appendix B Eqs. (B1)–(B5).

### Experimental/context source

K. Schneck et al. (SuperCDMS Collaboration),  
**“Dark matter effective field theory scattering in direct detection experiments,”**  
*Physical Review D* **91**, 092004 (2015), arXiv:1503.03379v3.

Schneck et al. use the same half-sum/half-difference isospin convention and the same 14-operator NREFT basis (excluding \(\mathcal O_2\)), following Anand et al. This paper is useful for experimental interpretation and is the paper originally supplied for this project.

### Internal provenance audit

`TARGETS_NUCLEAR_RESPONSES.md` records the project audit of nuclear-response sources, implementation references, and unresolved coefficient provenance. It does **not** override the analytic conventions frozen here.

---

## 2. Scope of this frozen convention

This convention describes **elastic, nonrelativistic WIMP–nucleus scattering** in the one-body Galilean-invariant EFT of Anand–Fitzpatrick–Haxton.

The formalism assumes:

- elastic scattering;
- the WIMP–nucleus interaction is obtained by summing one-body WIMP–nucleon interactions over nucleons;
- nuclear ground states are treated as having good parity and CP for the response reduction;
- the standard heavy-mediator/contact NREFT operator basis is used;
- the nuclear response functions are external nuclear-structure inputs.

This convention does **not** automatically cover:

- two-body/chiral nuclear currents;
- axion or dark-photon absorption;
- electron scattering;
- inelastic nuclear scattering;
- light-mediator propagators unless an explicit matching/propagator factor is added;
- alternative NREFT normalization conventions.

Those require explicit adapters or separate process models.

---

## 3. General units and software convention

The theory layer uses natural units,

\[
\hbar=c=1.
\]

Masses and momenta are expressed internally in GeV unless a function name explicitly states another unit. Velocities are dimensionless fractions of \(c\).

Project rule:

- target nuclear mass: `m_target_GeV`;
- WIMP mass: `m_chi_GeV`;
- recoil energy at public boundaries: normally `E_nr_keV`;
- momentum-squared at the nuclear-response boundary: `q2_GeV2`;
- Wilson coefficients: GeV\(^{-2}\).

Unit conversion is explicit at module boundaries. No implicit use of \(A\times0.9315\) GeV is part of the NREFT convention.

---

## 4. Momentum-transfer convention

Following Anand et al., the signed three-momentum transfer is

\[
\mathbf q=\mathbf p'-\mathbf p=\mathbf k-\mathbf k',
\]

where \(\mathbf p,\mathbf p'\) are incoming/outgoing WIMP momenta and \(\mathbf k,\mathbf k'\) are incoming/outgoing nucleon momenta [Anand et al., discussion surrounding Eq. (3) and Eq. (10)].

Thus the nuclear recoil momentum is \(-\mathbf q\). Scalar observables depend on \(q^2\), but the sign convention matters when comparing operator definitions across formalisms.

The NREFT operators use the dimensionless combination

\[
\frac{\mathbf q}{m_N},
\]

where \(m_N\) is the **nucleon reference mass used in the EFT normalization**. It is not the isotope mass and must never be replaced by \(m_T\).

The symbolic definition \(Q=q^2/m_N^2\) is frozen. The exact numerical value of the nucleon reference mass is deferred to the future kinematics/rate assembly. The particle-response kernel accepts \(Q\) directly and must not introduce a numerical \(m_N\).

The future rate layer must provide `m_N_reference_GeV` explicitly and record it for reproducibility; no hidden default is allowed. Changing this numerical value does not by itself require a new analytic `convention_id`, provided the mathematical definition \(Q=q^2/m_N^2\) and the Wilson-coefficient convention are unchanged. It does change numerical predictions and must therefore be explicit in analysis configuration/provenance and trigger the relevant numerical regression tests.

---

## 5. Transverse velocity

For WIMP–nucleon scattering,

\[
\mathbf v=\mathbf v_{\chi,\mathrm{in}}-\mathbf v_{N,\mathrm{in}},
\]

and Anand et al. define

\[
\mathbf v^\perp=\mathbf v+\frac{\mathbf q}{2\mu_N},
\]

with

\[
\mathbf v^\perp\cdot\mathbf q=0
\]

by elastic energy conservation [Anand et al., Eqs. (9)–(10)].

At nuclear level, the particle response functions depend on the WIMP–target center-of-mass transverse velocity \(\mathbf v_T^\perp\). For an elastic target of mass \(m_T\),

\[
\mu_T=\frac{m_\chi m_T}{m_\chi+m_T},
\]

and, with the momentum convention above,

\[
v_T^{\perp\,2}=v^2-\frac{q^2}{4\mu_T^2}.
\]

At fixed recoil energy,

\[
q^2=2m_T E_R,
\qquad
v_{\min}=\frac{q}{2\mu_T}.
\]

The package must not treat \(v_T^{\perp\,2}\) and \(v^2\) as interchangeable.

---

## 6. NREFT operator basis

The canonical operator basis is the Anand/Schneck basis

\[
\{\mathcal O_1,\mathcal O_3,\mathcal O_4,\ldots,\mathcal O_{15}\},
\]

with \(\mathcal O_2\) excluded.

The operators are

\[
\begin{aligned}
\mathcal O_1 &= 1_\chi 1_N,\\
\mathcal O_3 &= i\mathbf S_N\cdot\left(\frac{\mathbf q}{m_N}\times\mathbf v^\perp\right),\\
\mathcal O_4 &= \mathbf S_\chi\cdot\mathbf S_N,\\
\mathcal O_5 &= i\mathbf S_\chi\cdot\left(\frac{\mathbf q}{m_N}\times\mathbf v^\perp\right),\\
\mathcal O_6 &=\left(\mathbf S_\chi\cdot\frac{\mathbf q}{m_N}\right)\left(\mathbf S_N\cdot\frac{\mathbf q}{m_N}\right),\\
\mathcal O_7 &= \mathbf S_N\cdot\mathbf v^\perp,\\
\mathcal O_8 &= \mathbf S_\chi\cdot\mathbf v^\perp,\\
\mathcal O_9 &= i\mathbf S_\chi\cdot\left(\mathbf S_N\times\frac{\mathbf q}{m_N}\right),\\
\mathcal O_{10} &= i\mathbf S_N\cdot\frac{\mathbf q}{m_N},\\
\mathcal O_{11} &= i\mathbf S_\chi\cdot\frac{\mathbf q}{m_N},\\
\mathcal O_{12} &= \mathbf S_\chi\cdot(\mathbf S_N\times\mathbf v^\perp),\\
\mathcal O_{13} &= i(\mathbf S_\chi\cdot\mathbf v^\perp)\left(\mathbf S_N\cdot\frac{\mathbf q}{m_N}\right),\\
\mathcal O_{14} &= i\left(\mathbf S_\chi\cdot\frac{\mathbf q}{m_N}\right)(\mathbf S_N\cdot\mathbf v^\perp),\\
\mathcal O_{15} &= -\left(\mathbf S_\chi\cdot\frac{\mathbf q}{m_N}\right)\left[(\mathbf S_N\times\mathbf v^\perp)\cdot\frac{\mathbf q}{m_N}\right].
\end{aligned}
\]

Source: Anand et al., Eqs. (12)–(14); Schneck et al., Eq. (2).

\(\mathcal O_{16}\) is not independent; Anand et al. give

\[
\mathcal O_{16}=\mathcal O_{15}+\frac{q^2}{m_N^2}\mathcal O_{12}.
\]

---

## 7. Isospin convention — FROZEN

The canonical internal Wilson basis is

\[
c_i^0,\qquad c_i^1,
\]

with

\[
\boxed{c_i^0=\frac{c_i^p+c_i^n}{2},\qquad c_i^1=\frac{c_i^p-c_i^n}{2}}
\]

and therefore

\[
\boxed{c_i^p=c_i^0+c_i^1,\qquad c_i^n=c_i^0-c_i^1.}
\]

The isospin operators are

\[
t^0=1,\qquad t^1=\tau_3,
\]

with

\[
\tau_3|p\rangle=+|p\rangle,\qquad \tau_3|n\rangle=-|n\rangle.
\]

Thus the interaction is

\[
\mathcal H_\mathrm{NREFT}=\sum_{\tau=0,1}\sum_i c_i^\tau\,\mathcal O_i\,t^\tau.
\]

Source: Anand et al., Eqs. (15)–(18); Schneck et al., Sec. I.

### Important consequence

Pure proton coupling:

\[
c_i^n=0\quad\Longleftrightarrow\quad c_i^0=c_i^1=\frac{c_i^p}{2}.
\]

Pure neutron coupling:

\[
c_i^p=0\quad\Longleftrightarrow\quad c_i^0=-c_i^1=\frac{c_i^n}{2}.
\]

Equal proton/neutron coupling:

\[
c_i^p=c_i^n=c\quad\Longleftrightarrow\quad c_i^0=c,\quad c_i^1=0.
\]

No alternative full-sum convention is allowed implicitly. Any imported dataset or external code using a different convention requires an explicit adapter.

---

## 8. Wilson-coefficient dimensions

The nonrelativistic coefficients have dimensions

\[
[c_i^\tau]=E^{-2}.
\]

The package stores them canonically in \(\mathrm{GeV}^{-2}\).

Anand et al. introduce a weak scale

\[
m_v=(\sqrt{2}\,G_F)^{-1/2}\approx246.22\ \mathrm{GeV}.
\]

and sometimes write

\[
c_i=\frac{\tilde c_i}{m_v^2}
\]

to quote dimensionless coefficients [Anand et al., Eqs. (19)–(20)].

**Source typo:** Anand et al., PRC 89, 065501, Eq. (20) prints \((2G_F)^{-1/2}=246.2\ \mathrm{GeV}\), which is algebraically inconsistent. The numerical value corresponds to the standard \((\sqrt{2}\,G_F)^{-1/2}\) definition. This source typo does not affect our dimensionful Wilson-coefficient implementation.

**Project rule:** \(\tilde c_i\) is a presentation/input convenience only. The internal NREFT engine uses the dimensionful \(c_i^\tau\) in GeV\(^{-2}\).

No factor of \(m_v^2\) is hidden inside \(R_k\), \(W_k\), or the rate engine.

---

## 9. WIMP spin

The WIMP spin \(j_\chi\) is an explicit model input.

The response formalism is not frozen to \(j_\chi=1/2\). Terms involving WIMP spin carry factors of

\[
j_\chi(j_\chi+1).
\]

There must be no hidden global/default WIMP spin in the generic interaction kernel.

---

## 10. Particle/nuclear factorization

Define the spin-averaged nonrelativistic transition probability

\[
P_\mathrm{tot}\equiv\frac{1}{2j_\chi+1}\frac{1}{2J+1}\sum_{\mathrm{spins}}|\mathcal M|^2.
\]

The central Anand factorization is

\[
P_\mathrm{tot}=\frac{4\pi}{2J+1}\sum_{\tau,\tau'=0,1}\left[R_M^{\tau\tau'}W_M^{\tau\tau'}+R_{\Sigma''}^{\tau\tau'}W_{\Sigma''}^{\tau\tau'}+R_{\Sigma'}^{\tau\tau'}W_{\Sigma'}^{\tau\tau'}+\frac{q^2}{m_N^2}\,\mathcal C^{\tau\tau'}\right],
\]

where

\[
\begin{aligned}
\mathcal C^{\tau\tau'}={}&R_{\Phi''}^{\tau\tau'}W_{\Phi''}^{\tau\tau'}+R_{\Phi''M}^{\tau\tau'}W_{\Phi''M}^{\tau\tau'}\\
&+R_{\tilde\Phi'}^{\tau\tau'}W_{\tilde\Phi'}^{\tau\tau'}+R_{\Delta}^{\tau\tau'}W_{\Delta}^{\tau\tau'}\\
&+R_{\Delta\Sigma'}^{\tau\tau'}W_{\Delta\Sigma'}^{\tau\tau'}.
\end{aligned}
\]

Source: Anand et al., Eq. (40).

### Critical normalization rule

The explicit factor

\[
\frac{q^2}{m_N^2}
\]

multiplying the five composite/interference channels in Eq. (40) is part of the **contraction convention**. It must not be silently absorbed into the stored \(W_k\) or duplicated inside the rate engine.

The three channels without this external factor are

\[
M,\qquad \Sigma'',\qquad \Sigma'.
\]

The five channels with it are

\[
\Phi'',\qquad \Phi''M,\qquad \tilde\Phi',\qquad \Delta,\qquad \Delta\Sigma'.
\]

---

## 11. Canonical WIMP response functions \(R_k^{\tau\tau'}\)

For compactness define

\[
Q\equiv \frac{q^2}{m_N^2},\qquad V\equiv v_T^{\perp\,2},\qquad S_\chi\equiv j_\chi(j_\chi+1).
\]

Then the canonical response functions are a direct transcription of Anand et al., Eq. (38):

\[
R_M^{\tau\tau'}=c_1^\tau c_1^{\tau'}+\frac{S_\chi}{3}\left[QV\,c_5^\tau c_5^{\tau'}+V\,c_8^\tau c_8^{\tau'}+Q\,c_{11}^\tau c_{11}^{\tau'}\right],
\]

\[
R_{\Phi''}^{\tau\tau'}=\frac{Q}{4}c_3^\tau c_3^{\tau'}+\frac{S_\chi}{12}(c_{12}^\tau-Qc_{15}^\tau)(c_{12}^{\tau'}-Qc_{15}^{\tau'}),
\]

\[
R_{\Phi''M}^{\tau\tau'}=c_3^\tau c_1^{\tau'}+\frac{S_\chi}{3}(c_{12}^\tau-Qc_{15}^\tau)c_{11}^{\tau'},
\]

\[
R_{\tilde\Phi'}^{\tau\tau'}=\frac{S_\chi}{12}\left[c_{12}^\tau c_{12}^{\tau'}+Q\,c_{13}^\tau c_{13}^{\tau'}\right],
\]

\[
\begin{aligned}
R_{\Sigma''}^{\tau\tau'}={}&\frac{Q}{4}c_{10}^\tau c_{10}^{\tau'}\\
&+\frac{S_\chi}{12}\big[c_4^\tau c_4^{\tau'}+Q(c_4^\tau c_6^{\tau'}+c_6^\tau c_4^{\tau'})+Q^2c_6^\tau c_6^{\tau'}\\
&\qquad\qquad+Vc_{12}^\tau c_{12}^{\tau'}+QVc_{13}^\tau c_{13}^{\tau'}\big],
\end{aligned}
\]

\[
\begin{aligned}
R_{\Sigma'}^{\tau\tau'}={}&\frac18\left[QV\,c_3^\tau c_3^{\tau'}+V\,c_7^\tau c_7^{\tau'}\right]\\
&+\frac{S_\chi}{12}\bigg[c_4^\tau c_4^{\tau'}+Qc_9^\tau c_9^{\tau'}\\
&\qquad+\frac{V}{2}(c_{12}^\tau-Qc_{15}^\tau)(c_{12}^{\tau'}-Qc_{15}^{\tau'})+\frac{QV}{2}c_{14}^\tau c_{14}^{\tau'}\bigg],
\end{aligned}
\]

\[
R_{\Delta}^{\tau\tau'}=\frac{S_\chi}{3}\left[Qc_5^\tau c_5^{\tau'}+c_8^\tau c_8^{\tau'}\right],
\]

\[
R_{\Delta\Sigma'}^{\tau\tau'}=\frac{S_\chi}{3}\left[c_5^\tau c_4^{\tau'}-c_8^\tau c_9^{\tau'}\right].
\]

These functions belong to the **particle/interactions layer**. They must not contain isotope abundances, detector effects, halo-model parameters, or nuclear response coefficients.

---

## 12. Nuclear response functions \(W_k^{\tau\tau'}\)

The canonical nuclear responses are

\[
M,\quad \Sigma'',\quad \Sigma',\quad \Phi'',\quad \tilde\Phi',\quad \Delta,
\]

plus the ordered interference responses

\[
\Phi''M,\qquad \Delta\Sigma'.
\]

They are defined as bilinears of reduced nuclear matrix elements in Anand et al., Eqs. (39)–(41).

For a harmonic-oscillator representation, they are functions of

\[
y=\left(\frac{qb}{2}\right)^2.
\]

### Project storage rule

The package's public nuclear-response boundary accepts \(q^2\) in GeV\(^2\), and a dataset may internally use \(q^2\), \(q\), \(y\), or another explicitly declared source coordinate.

For an Anand-style dataset:

\[
y=\frac{q^2b^2}{4}.
\]

The oscillator parameter \(b\) is external, provenance-bearing nuclear input. There is no hidden \(b(A)\) default in production data.

### No hidden transformations

A stored \(W_k\):

- contains nuclear structure only;
- contains no Wilson coefficients;
- contains no isotope abundance;
- contains no halo factor;
- contains no detector response;
- does not automatically absorb the external \(q^2/m_N^2\) factors of Eq. (40).

Response-dataset metadata must carry

```text
convention_id = anand2014_prc89_065501_half_isospin_v1
```

before it can be contracted with the canonical \(R_k\) implementation.

---

## 13. Ordered isospin/interference entries

The contraction sums explicitly over

\[
\tau,\tau'\in\{0,1\}.
\]

The response layer must preserve ordered keys \((\tau,\tau')\). No implicit replacement

\[
(0,1)\leftrightarrow(1,0)
\]

is performed unless a source or explicit adapter establishes the required symmetry.

Likewise the mixed nuclear responses \(\Phi''M\) and \(\Delta\Sigma'\) are stored as the ordered objects defined in Anand et al., Eq. (41), not renamed or symmetrized ad hoc.

---

## 14. Differential cross section normalization

Using the nonrelativistic/Galilean transition probability \(P_\mathrm{tot}\) defined above,

\[
\boxed{\frac{d\sigma}{dE_R}=\frac{m_T}{2\pi v^2}\,P_\mathrm{tot}.}
\]

Source: Anand et al., Eq. (50), equivalently Appendix B Eq. (B4).

Here:

- \(m_T\) is the actual target-nucleus mass used by the target layer;
- \(v\) is the incoming WIMP speed in the laboratory frame;
- \(P_\mathrm{tot}\) uses the nonrelativistic normalization.

Do not insert the relativistic-amplitude conversion factor \((4m_\chi m_T)^2\) into the canonical rate path. Anand's Mathematica package offers that as an optional presentation convention; it is not our internal normalization.

---

## 15. Differential event rate

For target nuclei of a single isotope,

\[
\frac{dR}{dE_R}=N_T\,\frac{\rho_\chi}{m_\chi}\int_{v>v_{\min}}d^3v\,f_E(\mathbf v)\,v\,\frac{d\sigma}{dE_R},
\]

where \(N_T\) is the number of target nuclei per chosen detector-mass normalization.

Equivalently,

\[
\boxed{\frac{dR}{dE_R}=N_T\,\frac{\rho_\chi m_T}{2\pi m_\chi}\left\langle\frac{P_\mathrm{tot}(v,q)}{v}\right\rangle_{v>v_{\min}}.}
\]

Source: Anand et al., Eqs. (52)–(53) and Appendix B Eq. (B5).

### Mixtures

For a detector containing several isotopes, compute the rate for each isotope with its own \(m_T\), \(J\), nuclear-response dataset, and number of target nuclei per detector mass, then sum isotope rates.

No isotope abundance is folded into \(W_k\).

The existing target/composition layer is responsible for distinguishing number fractions from mass fractions and for explicit conversions between them.

---

## 16. Halo separation

The NREFT convention does not define a particular halo model.

The halo layer supplies the velocity averages required by the rate. Because different \(R_k\) contain different powers of \(v_T^{\perp\,2}\), the generic rate engine must support the required velocity moments rather than hard-code only the standard mean inverse speed \(\eta(v_{\min})\).

Halo parameters such as \(v_0\), \(v_\mathrm{esc}\), \(v_\mathrm{lab}\), and \(\rho_\chi\) are analysis assumptions, not NREFT conventions.

---

## 17. Mandatory \(\mathcal O_1\) normalization benchmark

Before implementing the full operator set in production, the generic engine must reproduce the standard SI result for \(\mathcal O_1\).

If only \(c_1^\tau\) are nonzero,

\[
R_M^{\tau\tau'}=c_1^\tau c_1^{\tau'},
\]

and every other \(R_k\) vanishes.

Therefore

\[
P_\mathrm{tot}^{(\mathcal O_1)}=\frac{4\pi}{2J+1}\sum_{\tau,\tau'}c_1^\tau c_1^{\tau'}W_M^{\tau\tau'}(y).
\]

At \(q\to0\),

\[
\sqrt{4\pi}\,M_{00;0}(0)=A,\qquad \sqrt{4\pi}\,M_{00;1}(0)=Z-N,
\]

so the coherent amplitude becomes

\[
c_1^0A+c_1^1(Z-N)=c_1^pZ+c_1^nN.
\]

Thus

\[
\boxed{P_\mathrm{tot}^{(\mathcal O_1)}(q=0)=(c_1^pZ+c_1^nN)^2.}
\]

For equal proton/neutron coupling,

\[
c_1^p=c_1^n=c\quad\Rightarrow\quad P_\mathrm{tot}(0)=c^2A^2.
\]

For a free proton target,

\[
\boxed{\sigma_p^{(\mathcal O_1)}=\frac{\mu_{\chi p}^2}{\pi}(c_1^p)^2.}
\]

Consequently,

\[
\sigma_T^{(\mathcal O_1)}(q=0)=\frac{\mu_{\chi T}^2}{\pi}(c_1^pZ+c_1^nN)^2.
\]

This is the primary absolute-normalization regression test for the future NREFT engine.

Required numerical validation:

\[
\text{dedicated SI/O}_1\;\approx\;\text{generic NREFT/O}_1\;\approx\;\text{independent external implementation},
\]

including absolute normalization, not merely spectral shape.

---

## 18. Software-layer contract

The intended dependency structure is

```text
Wilson coefficients + WIMP spin
            |
            v
    particle responses R_k
            |
            |        external nuclear data
            |                |
            |                v
            +---------> nuclear W_k
                         |
halo moments ---------->|
                         v
                    rate engine
                         |
                         v
                 truth recoil spectrum
                         |
                         v
                 detector response
                         |
                         v
                 statistics/sensitivity
```

### `interactions/`

Owns \(c_i^\tau\), WIMP spin \(j_\chi\), operator basis, \(R_k^{\tau\tau'}\), and the explicit convention ID. It does not own nuclear data, halo assumptions, detector effects, or backgrounds.

### `targets/`

Owns isotope identity/mass/spin, composition, external \(W_k^{\tau\tau'}\) datasets, source coordinates and provenance, and nuclear-response convention ID. It does not own Wilson coefficients or halo integrals.

### `rates/`

Owns compatibility checking of convention IDs, the \(R_kW_k\) contraction, explicit Eq. (40) \(q^2/m_N^2\) factors, differential-cross-section normalization, isotope summation, and halo-moment contraction. It must not reinterpret or silently convert conventions.

---

## 19. Known isolated issue: Mathematica `SetCoeffsNonrel`

The analytic equations in the peer-reviewed Anand et al. paper unambiguously define

\[
c_i^0=\frac{c_i^p+c_i^n}{2},\qquad c_i^1=\frac{c_i^p-c_i^n}{2}.
\]

However, Appendix B states that

```text
SetCoeffsNonrel[4,12.3,0]
```

produces

\[
c_4^p=c_4^n=6.15,
\]

which is inconsistent with the analytic inverse relation

\[
c_i^p=c_i^0+c_i^1,\qquad c_i^n=c_i^0-c_i^1
\]

if the input `12.3` is literally \(c_4^0\).

**Project decision:** this discrepancy does not alter the canonical analytic convention. The equations in the peer-reviewed paper and the matching Schneck convention are authoritative for WIMP-LEGEND.

Before using the historical Mathematica script as an absolute numerical validation target, inspect the versioned source to determine whether its input routine uses a rescaled user-facing coefficient or whether the documentation example is erroneous.

This is a validation-adapter issue, not a blocker for the theoretical implementation.

---

## 20. Explicitly not frozen by this document

The following are deliberately outside this theoretical convention and must be sourced/configured separately:

- exact numerical value of the nucleon reference mass \(m_N\);
- exact isotope nuclear masses;
- Ge isotope fractions for a specific LEGEND configuration;
- harmonic-oscillator parameters \(b\);
- numerical \(W_k\) polynomial/table coefficients;
- trusted validity ranges of nuclear-response datasets;
- local dark-matter density;
- Standard Halo Model parameter values;
- detector quenching, resolution, efficiency, threshold;
- background model;
- statistical method;
- UV-to-NREFT matching for a particular dark-matter model;
- two-body/chiral-EFT response extensions.

Changing one of these does not by itself change the NREFT convention ID unless it changes the mathematical definitions of \(c_i^\tau\), \(R_k\), \(W_k\), or their contraction.

---

## 21. Acceptance criteria for implementation

An implementation claiming compatibility with

```text
anand2014_prc89_065501_half_isospin_v1
```

must satisfy all of the following:

1. Use the half-sum/half-difference isospin convention exactly.
2. Use the operator basis and momentum sign defined above.
3. Store \(c_i^\tau\) internally in GeV\(^{-2}\).
4. Exclude \(\mathcal O_2\) from the canonical 14-operator basis.
5. Implement Eq. (38) literally for the eight particle-response combinations.
6. Implement Eq. (40) with the external \(q^2/m_N^2\) factors in the correct five channels.
7. Reject incompatible nuclear-response `convention_id` values.
8. Use the nonrelativistic transition-probability normalization of Eq. (B3).
9. Use \(d\sigma/dE_R=m_TP_\mathrm{tot}/(2\pi v^2)\).
10. Reproduce the \(\mathcal O_1\) coherent \(A^2\) limit.
11. Reproduce the dedicated SI implementation in absolute rate.
12. Independently cross-check at least one benchmark against a second implementation after explicitly mapping its conventions.

---

## 22. References

1. N. Anand, A. L. Fitzpatrick, W. C. Haxton,  
   *Weakly interacting massive particle-nucleus elastic scattering response*,  
   **Phys. Rev. C 89, 065501 (2014)**.  
   DOI: https://doi.org/10.1103/PhysRevC.89.065501

2. K. Schneck et al. (SuperCDMS Collaboration),  
   *Dark matter effective field theory scattering in direct detection experiments*,  
   **Phys. Rev. D 91, 092004 (2015)**, arXiv:1503.03379v3.  
   https://arxiv.org/abs/1503.03379

3. A. L. Fitzpatrick, W. C. Haxton, E. Katz, N. Lubbers, Y. Xu,  
   *The Effective Field Theory of Dark Matter Direct Detection*,  
   **JCAP 02 (2013) 004**, arXiv:1203.3542.

---

## 23. Change-control rule

This file is normative for the WIMP-LEGEND NREFT implementation.

If a future change modifies any of the following,

- isospin normalization;
- operator definition/sign;
- Wilson-coefficient dimensions;
- \(R_k\) definitions;
- \(W_k\) normalization;
- location of the explicit \(q^2/m_N^2\) factors;
- spin-average normalization;
- differential-cross-section normalization,

then the convention ID must be changed and an explicit migration/adapter must be provided.

External data updates that preserve these definitions do not require a new convention ID.
