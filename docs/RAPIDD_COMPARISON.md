# RAPIDD cross-check: first-stage results and tester handoff

## All-operator extension

Open [the verification notebook](../notebooks/rapidd_verification.ipynb) to
plot the saved comparison or rerun it. The expanded scan covers all 14
single-operator kernels in the frozen basis (O1 and O3–O15). O2 is outside
the basis. The original three cases match the poster workbench scope.
There are 840 kernel comparisons on the original grid and 640 matched-y
nuclear comparisons including the candidate mapping below.

With shared nuclear inputs, every nonzero shifted kernel on this grid
agrees within 3.75e-13 relative. Both implementations return zero for 336
other cases. The raw assembly reaches about 42.6% relative difference
for O5 and O8. This checks diagonal fixed-speed kernels; mixed operators,
absolute halo-integrated rates, folding and inference remain open.

**Correction to the first-stage scope:** RAPIDD's Ge-73 `pn_Fitz` table
contains `Pp/Pp` entries. Its O12 and O13 kernels use `Pp/Pp + Pt/Pt`,
which the extension maps to the canonical Phi-tilde-prime channel for
diagnosis. The pinned Ge table has no `Pt/Pt` entries. This mapping gives
shared-input kernel agreement, while the matched-y nuclear-input comparison
shows RAPIDD/AFH ratios of about 5 at y=0 for the Ge-73 p/n entries. The
neutron–neutron ratio reaches about 6.59 at y=0.1. The tabulated RAPIDD
Ge-73 `Pp/Pp` prefactor is `2*pi`; its normalization and polynomial
coefficients warrant separate review. The extension uses the raw table
values, with the candidate mapping explicitly recorded in the output.

## Status

Two reproducible checks have been completed against unmodified, compiled RAPIDD
C code: 560 nuclear-input evaluations and 180 shared-input O1/O8/O11 kernel
comparisons. **This is not yet a complete differential-rate, detector-folding,
exclusion-limit or surrogate validation.** No production physics was changed.

RAPIDD is a surrogate framework: Cerdeño, Cheek, Reid and Schulz,
[JCAP 08 (2018) 011 / arXiv:1802.03174v2](https://arxiv.org/abs/1802.03174v2).
Its numerical C calculator and its trained polynomial surrogate must be tested
separately. Agreement between implementations using related nuclear inputs is
not an independent nuclear-structure calculation.

## Frozen inputs and provenance

- WIMP-LEGEND baseline: `86cf0bd3beeb138cdc1e00db87992923720ace20`.
- RAPIDD repository: <https://github.com/cheekyparticle/RAPIDD_for_DM>.
- RAPIDD revision actually tested: `c36b290abfc8f68dc4e53445bc18c185f4b466b0`.
  This is **not** the earlier historical reference `02035a29...`.
- Actual C source: `lib/source/effFormFact.c`, SHA-256
  `b4aaa7b191ea012c71d29f8a4e33fc03850f98f850d559dce94370322db210a0`.
- AFH dataset manifest SHA-256:
  `b17df65326c24a46a4effef6f2b5881304c9c0f03792417c2635e921a6a76f33`.
- AFH original supplement SHA-256:
  `943a403e216633840887bbfcbec55c56e0c5813b86335b49868acb89e74a49b9`.
- Canonical production convention remains
  `anand2014_prc89_065501_half_isospin_v1`.

The RAPIDD `pn_Fitz` table identifies rounded Fitzpatrick form factors and warns
against its high-recoil use. It is an implementation source, not a primary
nuclear calculation. Its source comments cite arXiv:1203.3542 and
arXiv:1308.6288; our AFH export retains the original supplement density closure.
No external coefficients are copied into trusted production data. The JSON
comparison outputs contain diagnostic evaluations, not replacement datasets.

## 1. Nuclear inputs

All five stable Ge isotopes, four ordered p/n entries, seven channels, and
`y = 0, 0.001, 0.01, 0.1` were compared. The first-stage artifact excluded
`Phi_tilde_prime`; see the correction above. The seven channels are
M, Sigma-prime, Sigma-double-prime, Phi-double-prime, Delta, and the two mixed
channels. Source mixed-channel ordering is retained before canonical
transposition. Missing-entry C zeros are not independently certified by this
comparison.

RAPIDD defines recoil-to-y internally using
`m_A = Z*0.938 + (A-Z)*0.940 GeV` and its `B2` table. Consequently the script
chooses RAPIDD recoil energies that give the **same y**, rather than silently
comparing different momentum coordinates at the same recoil energy.
The comparison uses our stored source p/n polynomials directly, before the
canonical W transformation. There is no accidental factor-four or spin factor.

Maximum absolute relative differences among nonzero entries on this grid:

| Channel | Maximum |
|---|---:|
| M | 0.409% |
| Sigma-prime | 3.259% |
| Sigma-double-prime | 3.890% |
| Phi-double-prime | 7.925% |
| Delta | 3.159% |
| Phi-double-prime–M | 3.829% |
| Delta–Sigma-prime | 1.899% |

These are individual entries, not total rates or limits; cancellation can amplify
relative errors. Rounded coefficients alone need not guarantee uniform relative
accuracy near zeros. At y=0, RAPIDD's isoscalar coherent M sum is exactly A²;
our unrescaled source is smaller. RAPIDD/source minus one is respectively
0.044524%, 0.006090%, 0.016263%, 0.037912%, 0.058654% for
A=70,72,73,74,76. Do not repair the documented AFH closure deficit.

## 2. Shared-input particle/nuclear kernels

At each comparison point the C p/n responses are transformed into temporary
canonical W fixtures:

`W = (2J+1)/(4*pi) S F S^T`, `S = [[1,1],[1,-1]]`.

Mixed source matrices are transposed before this conversion. This is an explicit
comparison-only adapter, not automatic production convention conversion.
The unused Phi-tilde-prime entry is a labelled diagnostic zero. These fixtures
are pointwise constants and must never be used as physical response functions.

Cases: five isotopes; ER=1,10,30 keV; mchi=40 GeV; jchi=1/2;
(cp,cn)=(1,1),(1,0),(1,-0.7),(0,1); speed=700 km/s;
explicit m_N=0.938 GeV to match this C implementation. Unit-amplitude
coefficients are used only to compare homogeneous kernels. The frozen physical
Wilson convention is unchanged.

| Operator | Largest relative difference using raw C velocity assembly |
|---|---:|
| O1 | 7.62e-14 |
| O11 | 2.98e-13 |
| O8 | 42.622% |

The O8 difference is isolated from nuclear data. The tested C assembly combines
`FormFact_v0 + v² FormFact_v2`; `difRateGen.c` weights those terms with raw
inverse-speed and speed moments. Its `halo.c` integrands contain powers of lab
speed, not the subtraction v²-vmin². Our elastic convention uses v_perp².
For even Ge (Delta=0), the relative kernel excess is exactly
`vmin²/(v²-vmin²)` on this grid. Replacing v² by v_perp² **diagnostically**
restores O8 agreement to 8.07e-14. This is evidence of an apparent missing
velocity subtraction in the tested assembly, requiring independent review and
upstream clarification before an O8 limit comparison can be interpreted.
It is not a comparison of completed halo-integrated rates or an upstream fix.

## Reproduce

From WIMP-LEGEND root, with Python/NumPy/SciPy, GCC and Git:

```sh
git clone https://github.com/cheekyparticle/RAPIDD_for_DM ../RAPIDD_for_DM
git -C ../RAPIDD_for_DM checkout c36b290abfc8f68dc4e53445bc18c185f4b466b0
python analyses/validation/rapidd/compare_responses.py ../RAPIDD_for_DM \
  --output analyses/validation/rapidd/response_comparison.json
PYTHONPATH=src python analyses/validation/rapidd/compare_kernels.py ../RAPIDD_for_DM \
  --output analyses/validation/rapidd/kernel_comparison.json
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python -m pytest -q -W error \
  tests/test_rapidd_comparison_artifacts.py
```

Both scripts require a clean checkout at the exact pin, compile the original C
source into temporary directories, and ignore upstream tracked `.o` binaries.
No Mathematica, GSL, CMake or RAPIDD surrogate training is required for these
first-stage checks. Executed with GCC 14.2.0; 5 artifact tests pass in 0.02 s.
The artifact tests audit recorded output; rerunning the scripts is the live
independent C check.

## Independent tester: priorities before publication

1. Reproduce both JSON outputs from the clean pin and inspect mixed p/n ordering,
   spin factors and half-isospin mapping independently.
2. Review the O8 velocity issue in `effFormFact.c`, `difRateGen.c` and `halo.c`;
   include the historical RAPIDD pin to distinguish inherited behavior from a
   recent implementation change. Do not tune our convention to force agreement.
3. Build the full C calculator (requires GSL) and compare isotope-level dR/dER
   for O1 and O11 first. Explicitly align masses, oscillator lengths, fixed lab
   speed, SHM truncation, density, jchi, units and coefficient normalization.
   RAPIDD's dimensionless coefficients use c_p,n * (246.2 GeV)^2; this factor
   cancels its rate denominator and is not our dimensionful coefficient API.
   Do not use its Helm-based `difrate_SI_standard` as the AFH comparison.
4. Use identical nuclear data on both sides for algorithm validation, then
   separately quantify the accepted AFH versus rounded-input discrepancy.
   Include opposite-sign couplings, zeros, kinematic boundaries and mixed ops.
5. Inspect natural-Ge fractions and per-kg mass weights; don't assume equal
   composition semantics. Then match response, acceptance and statistics before
   comparing limits. RAPIDD takes abs(total rate) in one path; signed raw
   contributions must also be inspected so errors are not hidden.
6. Validate the surrogate only within its trained mass/coupling/halo domain,
   against its own exact calculator first. Report interpolation errors separately.

Publication gates still open: full rate/unit normalization, detector-folded
counts, inference, surrogate errors, nuclear-input systematics and independent
review. No publication-level acceptance is claimed by this first-stage report.
