# GE-NUCLEAR-DATA-001: AFH germanium supplement datasets

Version: `afh2014_ge_obdme_v1`. Five ground-state isotopes: Ge-70, 72, 73,
74, 76. Each has all eight channels and all four **ordered** isospin entries.
This is the one-body harmonic-oscillator calculation distributed with Anand,
Fitzpatrick and Haxton, *Phys. Rev. C* **89**, 065501 (2014),
DOI [10.1103/PhysRevC.89.065501](https://doi.org/10.1103/PhysRevC.89.065501).
It is not an Elastic, JUN45, jj44b or chiral two-body dataset.

## Public API and explicit analysis inputs

```python
from wimp_legend.targets import (
    germanium_afh2014_momentum, load_germanium_afh2014,
)

# ge76 is an existing Isotope with externally sourced bare-nucleus mass/spin.
# Calling this helper explicitly selects the supplement's b(A) and unit inputs.
momentum = germanium_afh2014_momentum(ge76.A)
dataset = load_germanium_afh2014(ge76, momentum=momentum)
W01 = dataset.evaluate(channel="Phi_double_prime_M", tau=0, tau_prime=1,
                      q2_GeV2=[0.0, 0.001, 0.01])
```

The loader returns the existing immutable `NuclearResponseDataset`. Its isotope
is the caller's object, including mass provenance; Ge identity and ground-state
spin must agree with the supplement (J=0 for even A, J=9/2 for A=73,
`IsotopicSpins`, line 947). No nuclear masses, abundance or enrichment values
are introduced. It performs no halo, WIMP, composition, detector or rate
calculation. Scalar/vectorized evaluation and ordered missing-key behavior are
those of the existing interface.

`momentum` is mandatory and must be a sourced `MomentumVariable("y", ...)`.
The explicit helper reproduces `SetIsotope`, line 1479:

\[
b_{\rm fm}=\sqrt{\frac{41.467}{45 A^{-1/3}-25 A^{-2/3}}},\quad
b_{\rm GeV^{-1}}=5.0677 b_{\rm fm},\quad y=q^2b^2/4.
\]

The number 5.0677 is the supplement's rounded `Femtometer` definition (line 17),
not a replacement for the package's physical constants. The helper records
both values and their source. A different sourced b (including a different fm
conversion) can be explicitly passed via `MomentumVariable`; no polynomial is
refitted. At equal y all such mappings give equal W, while fixed-q predictions
can change. The printed Fitzpatrick discussion's b=2.108 fm is **not** substituted
for the supplement's default. The frozen analytic convention does not change
when b changes; analysis provenance and numerical regression results must.

## Source, extraction and precision

The supplied `dmformfactor-prc(1).m`, 94,031 bytes / 1,634 lines, is pinned to
SHA-256 `943a403e216633840887bbfcbec55c56e0c5813b86335b49868acb89e74a49b9`.
`pfGe70`, `pfGe72`, `pfGe73`, `pfGe74`, `pfGe76` are at lines 373–524,
with row counts 15, 15, 83, 15, 15. Each exported density row retains its
original decimal string, quantum numbers, index and source line. Records are
`{Psi,K,T,N_bra,2j_bra,N_ket,2j_ket}`, N the major oscillator number. The written
56Ni core is included once. The package does not solve a shell-model problem.

`tools/ge_nuclear_data/afh_symbolic.py` transcribes the source single-particle
operators and harmonic-oscillator radial integrals (lines 803–922). Terminating
hypergeometrics are expanded exactly; angular recouplings use exact SymPy
3j/6j/9j symbols. `regenerate.py` parses only the checksum-pinned Ge rows and
spins. All written decimal OBDMEs are exact rationals during reconstruction;
coefficients are exported as decimal strings with 40 significant digits and
read as float64 by the production loader. The loader marks these derived
coefficient records as `kind="implementation"`, retaining the pinned primary
supplement and original publication as upstream provenance. Forty exported digits do not improve
the precision of the original, truncated source decimals. No fitting, numerical
interpolation, `Chop`, particle-number repair or deletion of small coefficients
is used. Exact algebraic selection-rule zeros are stored as the polynomial `[0]`.

The GCN2850 model label is supported by the supplied audit and reproduction of
Fitzpatrick's published form factors. The label is not written in the supplement
itself. The source comment at line 1470 points to arXiv:1203.3542v2; checks here
use the corrected v3 and its published JCAP normalization. Original APS download
date/account were not supplied; attachment acquisition is recorded as
2026-10-19, without inventing a historical download record.

## Provenance classification of numerical inputs

| Quantity | First appearance in supplied material | Role and traceability |
|---|---|---|
| Ge OBDMEs, including core | `pfGe70` ... `pfGe76`, lines 373–524 | Primary journal supplement nuclear-structure input; linked to Fitzpatrick's original calculation by source comment and published-table reproduction. Original shell-model run files/hash are not supplied. |
| A, Z, spin | Density-table names, `IsotopicSpins`, line 947 | Source isotope/property data; the supplement does not supply a separate experimental spin reference. |
| 41.467, 45, 25 in b(A) | `SetIsotope`, line 1479; also Fitzpatrick Appendix A.3 | Source model/oscillator prescription, explicitly selected; not measured isotope masses or a validated momentum range. |
| 5.0677 GeV^-1/fm | `Femtometer`, line 17 | Source rounded unit conversion, preserved by the explicit helper. |
| Exported W polynomial coefficients | This reconstruction of source rows and operators | Derived numerical data, not a second nuclear calculation. Each entry records original response index/isospin ordering and the documented normalization map. |
| Isotope mass, abundance, enrichment | Not supplied as trusted nuclear input | Caller input only; no numerical value imported from source `SetM`, other packages or the audit. |

The supplied comparison document is an audit, not primary nuclear-structure
data. No numerical Elastic, JUN45, jj44b or alternative implementation input is
imported. Nuclear-response coefficients remain separate from quantities the
existing NREFT engine computes from Wilson coefficients and kinematics.

## Explicit convention conversion

Canonical convention: `anand2014_prc89_065501_half_isospin_v1`.
`docs/NREFT_CONVENTIONS.md` and all production rate formulas are unchanged.
With source p/n form factors F (the `FF` routine, lines 959–977),

\[
S=\begin{pmatrix}1&1\\1&-1\end{pmatrix},\qquad
W_{\rm canonical}(y)=\frac{2J+1}{4\pi}e^{-2y}S F(y)S^T.
\]

This follows directly from cp=c0+c1, cn=c0-c1. The source's full-sum
isospin output instead uses one quarter of this expression. Thus canonical W
is **four times** the source `ResponseNuclear` W, with the mixed transpositions
below. There is no runtime guessing or automatic convention conversion.

| Source index | Stored canonical channel | Source operator order | Isospin action |
|---|---|---|---|
| 1 | `M` | M, M | unchanged |
| 2 | `Sigma_double_prime` | Sigma'', Sigma'' | unchanged |
| 3 | `Sigma_prime` | Sigma', Sigma' | unchanged |
| 4 | `Phi_double_prime` | Phi'', Phi'' | unchanged |
| 5 | `Phi_tilde_prime` | Phi-tilde', Phi-tilde' | unchanged |
| 6 | `Delta` | Delta, Delta | unchanged |
| 7 | `Phi_double_prime_M` | M, Phi'' | transpose (tau,tau') |
| 8 | `Delta_Sigma_prime` | Sigma', Delta | transpose (tau,tau') |

The transpositions follow the explicit source assignments, lines 995–996, and
canonical Eq. (41). Source-order p/n polynomials are retained in the JSON for
audit. They must not be contracted as canonical W. Within each channel entries
are ordered 00, 01, 10, 11, including zeros. The loader verifies the exact 32-key
set and hashes before constructing the dataset.

Stored W is the dimensionless reduced-matrix-element bilinear. The nuclear
operator's spherical-harmonic normalization is retained. The outer
4*pi/(2J+1), the five external q²/m_N² factors, WIMP coefficients and any
interference multiplicities remain in the frozen contraction. No m_N is
introduced here.

## Particle-number closure

At zero momentum the canonical O1 result is
`(cp*Z_eff + cn*N_eff)**2`, where the effective particle numbers are those
reconstructed from the supplied OBDMEs. They are not silently changed to integer
Z and N. For equal cp=cn=1:

| A | Effective Z | Effective N | P_tot(0)/A² |
|---|---|---|---|
| 70 | 31.997784920972 | 37.986636768915 | 0.999554955 |
| 72 | 31.999724641172 | 39.998082924996 | 0.999939100 |
| 73 | 31.999380070899 | 40.994684647114 | 0.999837396 |
| 74 | 31.998707487887 | 41.987269010945 | 0.999621022 |
| 76 | 31.998117684705 | 43.979603703472 | 0.999413807 |

Full-precision values are in each JSON. The deficit is documented in the
supplied audit; its physical/numerical origin remains unresolved. It is a source
property, not a different canonical normalization. Tests check this residual,
all ordered coherent entries, proton-only/neutron-only and unequal couplings,
and the existing generic NREFT contraction without changing that contraction.

## Domain and verification scope

All 160 entries use `exp(-2y) P(y)` on the mathematical domain **y>=0**, with no
invented upper bound. Negative q² and source coordinates are rejected. There is
**no source-backed finite physically validated range**, so `validated_range`
is None. This is neither a claim that the one-body nuclear model is valid at
arbitrary momentum nor a claim of unlimited float64 numerical accuracy. Users
must state their analysis range; future sourced validity restrictions must be
recorded separately from the mathematical domain.

Verification supplied here:

- 465 individually printed coefficients, all five isotopes and every printed
  channel, agree with Appendix A.3 of Fitzpatrick et al., *JCAP* 02 (2013) 004,
  [arXiv:1203.3542v3](https://arxiv.org/abs/1203.3542v3), within half a unit of
  the last printed significant digit. Omitted higher powers are retained.
  Alternating rounded coefficients are not used as a high-precision finite-y
  function oracle; cancellation makes that comparison misleading.
- `verify_source.py` executes the **original** pinned single-particle expressions
  in Mathics3 10.0.1 at y=1/4 and y=1. All 270 distinct Ge operator/orbit/rank
  combinations are included (540 values). Assembly from those references checks
  all 32 responses at both y for each isotope (320 values). Largest scaled
  relative response discrepancy observed: about 4.2e-13.
- Mathics' source NineJ sum bound is shortened from 50 to 5. The source SixJ
  triangle `(1,l,g)` requires `g<=l+1<=5` for these orbits, so no nonzero term is
  lost. All source operator expressions remain unchanged. Mathics and the Python
  port share a SymPy angular backend. This is a separate execution/transcription
  check, not commercial Mathematica certification or independent shell-model
  validation. Phi-tilde' has no corresponding printed Fitzpatrick table;
  its check is source replay plus exact algebra/selection rules.
- Exact even-isotope zeros, completeness/order, p/n transformation, mixed
  asymmetry, self-response PSD, shapes/immutability, explicit oscillator/mass
  handling, unknown validity and tampered-file rejection have dedicated tests.

No entry is substituted from another nuclear calculation. No channel is left
missing or populated with a guessed coefficient. The export's machine-readable
`verification_status="symbolic_source_reconstruction"` describes reconstruction,
not a physical validity certificate; the reference fixtures and checks above
state its actual verification scope.

## Reproduce and verify

From the repository root, with the original uploaded supplement available:

```bash
python -m pip install '.[test,nuclear-data-regeneration]'
python tools/ge_nuclear_data/regenerate.py \
  --supplement /path/to/dmformfactor-prc.m \
  --output /tmp/ge-afh2014-regenerated
diff -ru src/wimp_legend/targets/data/ge_afh2014_v1 /tmp/ge-afh2014-regenerated
python tools/ge_nuclear_data/verify_source.py \
  --supplement /path/to/dmformfactor-prc.m \
  --output /tmp/ge-afh2014-source-replay.json
python -m pytest -q -W error
```

Use SymPy 1.14.0 and Mathics3 10.0.1, pinned in the optional development extra.
The source replay is machine arithmetic: different interpreter/platform versions
may differ in final bits, so compare reference values with the documented test
tolerances, not byte identity. The analytic dataset export is deterministic and
should match byte for byte, including the manifest, with the pinned version.
Production installation requires neither SymPy, Mathics nor Mathematica.

`manifest.json` records per-isotope hashes and generator hashes; the loader pins
its own manifest hash. Changing data requires a new reviewed version and hash
update. Supporting readme/audit hashes are recorded in every isotope artifact.
See `GE_NUCLEAR_DATA_HANDOFF.md` for the independent tester's checks and the
remaining source/model questions.
