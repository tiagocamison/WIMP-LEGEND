# GE-NUCLEAR-DATA-001 — independent scientific acceptance

Audit date: 2026-10-19. Production code and repository tests were not modified. Expected responses were generated from the supplied primary density matrices with a separate numerical implementation, before reviewing the Coder's extraction implementation. Coder fixtures were not used as numerical oracles.

| Area | Verdict | Scientific basis / condition |
| --- | --- | --- |
| Dataset provenance and extraction | **PASS WITH CONDITIONS** | Supplied primary bytes, all 143 density rows, source locations, dataset hashes and numerical extraction agree. Original APS acquisition history and underlying shell-model artifacts remain unavailable; GCN2850 is an inferred identification, not a literal source label. Source closure defects remain. |
| Convention and interference transformations | **PASS** | Independent half-sum derivation requires ×4; both mixed channels require ordered-index transposition. Asymmetric entries and physical contractions agree. |
| Individual response-channel correctness | **PASS WITH CONDITIONS** | All eight channels independently reproduced for all five isotopes, including every exact-zero entry. Conditions: source closure defects, nominal one-body model scope, unknown physical validity cutoff, and float64 errors at diffraction nodes. |
| Integration with the physical rate engine | **PASS WITH CONDITIONS** | Independent contraction, absolute units, speed integration and recoil benchmarks agree. A prior architectural test fails because the convention ID is duplicated; tiny negative node responses can propagate to rates. No new ordinary-recoil normalization defect found. |

Nominal AFH one-body poster curves are supportable under these conditions. This is not an unconditional certification of nuclear-model accuracy or of a completed detector sensitivity calculation. No response channel remains unverified as an extraction from the supplied calculation. The nuclear Hamiltonian/wavefunctions, model uncertainty and a source-backed physical momentum cutoff have not been independently validated.

## Exact source identity

| Item | Identity |
| --- | --- |
| Reported Coder commit | `b2e6a423cb02acc64611d4e52dec047787345b85` — not retrievable from the remote |
| Exact tested, available commit | `c0283c1def642db3500184b68fcb042a323392b0` |
| Tested tree | **`94c524dfd5739feb648329f1154402f1b4844e30`** |
| Delivery baseline | `93a8ed143578ad4b6a74bd636ada2943fa7e0ecb` |
| Attached diff SHA-256 | `5c517075dbab55e689af60dba77d9f5fc99e29dafba8cf4c65a9fbd4890db891` |
| Original supplement | `dmformfactor-prc(1).m`, 94,031 bytes, 1,634 lines |
| Supplement SHA-256 | **`943a403e216633840887bbfcbec55c56e0c5813b86335b49868acb89e74a49b9`** |
| Frozen conventions SHA-256 | `a2c57658add45b273ca082d8b01fbab7dd739f8625d6f04a438d619d4cf385cd` |

Applying the supplied diff to the stated baseline in a separate temporary Git index produced exactly the tested tree. Thus the delivered file content is certified independently of the unavailable Coder commit object. The fetched poster tip `fcad8a4f415e28e830b011abbe7145bfd67a9f99` differs only by `.gitattributes`; it was not substituted for the delivered tree. Final tested checkout was clean. Full dataset hashes are in `SOURCE_IDENTITY.json` in the evidence archive.

The installation readme checksum is `cc844dc64f5b0cd255d2e8a21e219281634ee1e7d064439211c5eeb99e145648`; the supplied comparison report checksum is `230582e778f9653c68a22cc10e7f57461f2d4ea6125ab846a86b8fd81a2ba18f`. Both match dataset provenance. This establishes identity with the supplied files, not an independently recovered APS download. The Elastic comparison report is supporting context, not a numerical oracle or proof of identity with Elastic files.

As a separate provenance check, the unmodified Coder regeneration script was executed with pinned SymPy 1.14.0 in an isolated export environment. All five JSON datasets and the manifest were byte-identical to delivery. That verifies reproducibility of extraction, not physical correctness; the independent tests below supply the separate physics/numerical evidence. Forty exported digits do not improve the physical precision of the original density decimals.

## Independent references and results

`ge_reference.py` parses the original Mathematica density records directly. Angular 3j/6j symbols use factorial Racah sums with exact rational summation; 9j symbols use contracted 6j sums. Neither SymPy nor Mathics is imported. Normalized harmonic-oscillator wavefunctions and their derivatives are integrated directly with generalized Gauss–Laguerre quadrature, replacing the supplement/exporter's hypergeometric polynomial construction. The primary operator equations define the six single-particle operators; independently assembled multipoles define the eight responses. No production coefficients generate expected values.

Checks include 96-versus-192-node convergence, five separate 60-decimal mpmath radial integrations with directly differentiated wavefunctions, analytic 0s integrals, scalar charge and spin-orbit traces, spin longitudinal/transverse limits, and conserved total angular momentum. The model equations and source densities are shared physical inputs; this is not a new many-body nuclear calculation.

All 160 ordered canonical entries were compared at ten values of y: 0, 0.0001, 0.01, 0.1, 0.25, 0.5, 0.75, 1, 1.5, 2: **1,600 numerical comparisons**. The exported p/n polynomials were separately compared against the independent p/n reconstruction at the same points. Maximum 96/192-node absolute difference was `3.64e-12`, on responses of order thousands. Maximum production/reference channel-scaled error was `3.24e-13`. Scaling here divides by the largest absolute entry in that channel matrix at that isotope/y, with a `1e-10` floor; it does not imply relative accuracy at a zero.

| Channel | Scientific verdict | Nonzero isotopes | Maximum channel-scaled error | Additional physical check |
| --- | --- | --- | ---: | --- |
| M | PASS WITH CONDITIONS | All five | 3.24e-13 | Coherent charge trace; exact source deficit retained |
| Sigma double prime | PASS WITH CONDITIONS | 73 | 1.82e-15 | Spin limit; published neutron response |
| Sigma prime | PASS WITH CONDITIONS | 73 | 2.00e-14 | W(Sigma prime,0)=2 W(Sigma double prime,0) |
| Phi double prime | PASS WITH CONDITIONS | All five | 1.05e-13 | Independent L·sigma monopole trace; published response |
| Phi tilde prime | PASS WITH CONDITIONS | 73 | 2.41e-15 | Rank-2 q=0 relation; independent radial/angular evaluation |
| Delta | PASS WITH CONDITIONS | 73 | 8.49e-15 | Orbital-angular-momentum limit; published neutron response |
| Phi double prime M | PASS WITH CONDITIONS | All five | 8.78e-14 | Ordered mixed trace and published asymmetric p/n pair |
| Delta Sigma prime | PASS WITH CONDITIONS | 73 | 4.25e-14 | Ordered mixed response and published asymmetric p/n pair |

The four J=0 isotopes have only M, Phi double prime and their interference. The remaining five channels are exact zeros: **80 zero entries**, not numerical thresholding. Ge73 includes every density block K=0…9 for T=0,1. Independent parity/rank gates give even K for M/Phi, even K≥2 for Phi tilde prime, and odd K for both spin responses and Delta.

At q=0, the Ge73 Phi response includes a rank-2 term: the scalar spin-orbit trace alone is insufficient. Independently, `W_Phi(0) - (2/3) W_PhiTilde(0)` equals the outer product of the scalar spin-orbit traces. This relation passes. No second published numerical table for Phi tilde prime was found; its finite-q certification is against the supplied primary operator calculation, with an independent numerical backend.

Fitzpatrick et al., [arXiv:1203.3542v3, Appendix A.3](https://arxiv.org/html/1203.3542v3), supplies an additional external check. Eight manually transcribed Ge73 polynomials (62 rounded coefficients) were checked at y=0,0.01,0.25 against conservative two-significant-figure rounding envelopes. All 24 comparisons fit; the largest discrepancy uses 96.67% of its envelope. These printed functions are not precision numerical references: e.g. neutron Sigma prime at y=0.25 is `0.07128472690` independently, versus `0.08301035008` from printed coefficients, within an envelope of `0.03445572858`. Their limited precision must not be mistaken for an extraction discrepancy. The paper's b=2.108 fm is distinct from the supplement's explicitly selected A-dependent b(A); comparisons were made in y.

## Conventions and mixed-channel ordering

Let `S=[[1,1],[1,-1]]`. Frozen coefficients obey `c_half=(cp+cn,cp-cn)/2`, whereas the supplement's isospin setter stores `c_full=2 c_half`. Invariance of the bilinear implies `W_half=4 W_full`. Equivalently, canonical W is `(2J+1)/(4pi) S F_pn S^T`, including the source `exp(-2y)` factor when F is polynomial-only.

Source channel 7 is M first, Phi second; channel 8 is Sigma prime first, Delta second. Its first isospin index belongs to its first operator. Therefore canonical `W_PhiM[t,u]=4 W_source,MPhi[u,t]` and `W_DeltaSigma[t,u]=4 W_source,SigmaDelta[u,t]`. There is no extra interference factor and no symmetrization.

For Ge73 at y=0.25, independent canonical entries are:

| Channel | (0,1) | (1,0) |
| --- | ---: | ---: |
| Phi double prime M | 17.8570675913004 | -92.7420985162104 |
| Delta Sigma prime | 0.186005030230544 | 0.181073213396692 |

Production agrees, including signs and asymmetry. Arbitrary unequal left/right isospin vectors and interference-sensitive rate contractions also agree.

## Physical input defects: retain, document, do not renormalize

The coherent O1 identity is `P(0)=(cp Zeff+cn Neff)^2`, not the ideal integer-charge expression. Independent closure gives:

| A | Zeff | Neff | P(0)/A² for cp=cn=1 |
| --- | ---: | ---: | ---: |
| 70 | 31.9977849210 | 37.9866367689 | 0.999554954953 |
| 72 | 31.9997246412 | 39.9980829250 | 0.999939099987 |
| 73 | 31.9993800709 | 40.9946846471 | 0.999837396145 |
| 74 | 31.9987074879 | 41.9872690109 | 0.999621022368 |
| 76 | 31.9981176847 | 43.9796037035 | 0.999413806672 |

Equal-coupling SI rate deficits reach **0.05862%**. This small bound does not apply near isospin cancellation. For Ge76 the q=0 cancellation ratio shifts from ideal `cn/cp=-0.7272727273` to `-0.7275672128`; physical interpretation near such blind spots needs particular care. The supplied comparison report's statement that each neutron number closes to better than 3e-4 is too strong: Ge76 has a `4.64e-4` relative neutron deficit. The underlying deficit origin is unresolved.

A further independent source check finds `Jz=4.49904569645` for Ge73 instead of 4.5. Reconstructing total J from q=0 spin/orbital responses agrees with the direct diagonal K=1 density trace: reduced J `15.72879644921`, ideal `15.73213272255`, relative deficit **2.12067455e-4**. A half-last-decimal error bound on these diagonal densities is only `6.28e-7` in reduced J, versus the actual `0.00333627` deficit. Ordinary rounding of the written K=1 numbers cannot explain it. The strict conserved-J test fails and is retained. This is a defect/limitation of the supplied density input, not a newly introduced importer normalization error. No density or response was rescaled to hide either deficit.

## Integration, units and numerical limitations

Forty single-isotope contraction/rate cases cover O1, O4, O1+O3, O4+O5, O8+O9, O11+O12, O13 and a dense fourteen-operator configuration. They use independently assembled frozen particle bilinears and an analytic speed distribution, explicit synthetic target masses, `q²=2mT E_keV 1e-6`, and `Q=q²/m_N_reference²`. The nuclear spin factor is applied once; the five composite/interference channels receive one external Q. Absolute cross sections and the GeV-to-cm, seconds-to-days and GeV-to-keV factors agree. Kinematically inaccessible recoils are zero.

Eighteen additional halo benchmarks use v0=220, ve=232, vesc=544 km/s, mchi=50 GeV and rho=0.3 GeV/cm³. A separate 60-decimal speed-moment integration supplies expectations. Maximum rate discrepancy was `8.89e-16` relative. These exercise source M, spin and spin-orbit responses at physical recoil energies 1 and 10 keV. All coefficients are explicitly artificial benchmark choices; these are not sensitivity predictions.

For illustration, at 10 keV with `c_op^0=1e-4*(1+op/20)` and `c_op^1=-4e-5*(1-op/20)` GeV^-2, the independent Ge73 rates are `22096.85841149` (O1), `0.511546703987` (O4), and `15.31984022220` (O11+O12), in events/(kg day keV), with identical production results at displayed precision. All parameters are in the executable tests and `physical_rates.json`.

Benchmark masses are explicitly approximate: [NIST atomic masses](https://pml.nist.gov/cgi-bin/Compositions/stand_alone.pl?ascii=ascii) converted using [CODATA 2022](https://physics.nist.gov/cuu/Constants/Table/allascii.txt), `mT≈m_atom-32me`, omitting positive electronic binding energy/c². That approximation is recorded, not promoted to a precision bare-nucleus mass input. The loader correctly retains caller-supplied masses. Actual poster composition, masses, halo and detector inputs are outside this dataset acceptance.

Explicit alternate oscillator lengths map the same q² to the correct y. Negative/nonfinite q², wrong isotope/spin or momentum mapping and overflow-scale evaluation are rejected. Mathematical support is y≥0 without a finite upper boundary; absence of a validated physical range is exposed in metadata. Testing through y=2 establishes numerical agreement, not nuclear-model validity through that range.

Two strict numerical-boundary cases fail: expanded float64 M polynomials become negative at Ge70/72 diffraction zeros. At y=`1.078609041158255`, Ge70 W00 is `-6.57373634e-15` versus an independent squared amplitude of approximately `1.72e-30`; Ge72 gives `-1.88097143e-14` at y=`1.102269641831610`. These are tiny cancellation errors, not wrong physics coefficients, but a rate can inherit the negative sign. Relative accuracy at nodes or very near isospin blind spots is not certified. The positive residuals at Ge74/76 nodes also illustrate loss of relative precision. No clipping, test relaxation or production repair was applied.

## Regressions and reproducibility

Environment: Python 3.12.14; NumPy 2.3.5; SciPy 1.17.0; pytest 9.1.1; mpmath 1.4.1. The independent evaluator does not require the Coder's optional SymPy/Mathics stack.

| Suite | Result |
| --- | --- |
| Unchanged repository suite, warnings as errors | **830 passed**, 83.23 s |
| Unchanged external pre-rate, rate, generic analysis and background/SI audit suites | **652 passed, 1 failed**, 194.81 s |
| New independent scientific suite | **571 passed, 1 failed**, 11.57 s — failure is conserved-J closure of the primary input |
| Separate strict diffraction-node nonnegativity suite | **2 passed, 2 failed**, 0.49 s |

The external failure is `test_dependency_boundaries`: production now contains two literal copies of `anand2014_prc89_065501_half_isospin_v1`, in `interactions/wilson.py` and `targets/germanium.py`, where the earlier audit requires one authoritative definition. Their strings currently match; no numerical drift is demonstrated. This is a minor architectural implementation regression, not a nuclear-physics disagreement. The unchanged tests were neither edited nor excluded. The detector-remediation boundary suite was not rerun: detector algorithms are unchanged, and its previously documented remote Gaussian-tail limitation remains outside this milestone.

The evidence archive contains runnable tests, primary bytes, the delivered diff/context, unchanged external tests and their SHA-256 records, exact source identity, full numerical CSV/JSON comparisons, and raw logs. Strict failing checks remain failing. The README gives commands and expected exit statuses; no production changes are needed to reproduce this verdict.
