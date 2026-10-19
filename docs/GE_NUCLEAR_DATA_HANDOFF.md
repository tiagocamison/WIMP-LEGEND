# GE-NUCLEAR-DATA-001: independent verification handoff

Branch: `poster`; baseline commit `93a8ed143578ad4b6a74bd636ada2943fa7e0ecb`.
The delivery `SOURCE_IDENTITY.json` supplies the exact resulting commit and tree
without a self-referential hash inside the commit.

## Public boundary

`wimp_legend.targets.germanium_afh2014_momentum(A)` explicitly selects the source
oscillator prescription. `load_germanium_afh2014(isotope, *, momentum)` returns
an existing `NuclearResponseDataset`, with all 32 ordered keys. Supply sourced
bare-nucleus masses externally. No composition assumptions enter W.

Preserved analytic convention:
`anand2014_prc89_065501_half_isospin_v1`.
Normative document SHA-256 (unchanged):
`a2c57658add45b273ca082d8b01fbab7dd739f8625d6f04a438d619d4cf385cd`.

## Expected mathematical properties

1. Source F in p/n: `W=(2J+1)/(4*pi)*exp(-2y)*S@F@S.T`,
   `S=[[1,1],[1,-1]]`. This is four times the full-sum source W. Transpose the
   ordered isospin matrix of source channels 7 and 8 to canonical Phi''M and
   DeltaSigma'. Test unequal isospin entries; symmetrizing would hide this bug.
2. `4*pi/(2J+1) * c.T @ W_M(0) @ c` equals
   `(cp*Z_eff + cn*N_eff)**2`, with half-sum c and **source effective counts**.
   Exact integer A² is deliberately not enforced. Closure residuals are in
   `GE_NUCLEAR_DATA.md` and the exported normalization records.
3. Even Ge isotopes have exactly zero Sigma', Sigma'', Phi-tilde', Delta,
   DeltaSigma' entries, including all four isospin orders. Ge-73 has all eight
   channels. Self-channel matrices are symmetric PSD; mixed matrices need not be
   symmetric or positive. No external Q or Wilson factor belongs in these data.
4. All responses are `exp(-2y)*P(y)` with ascending coefficients, source written
   decimals treated exactly during extraction, 40 significant digits exported,
   float64 at runtime. Tiny coefficients must not be discarded merely because
   the printed appendix omitted them.
5. Mathematical domain: finite nonnegative y, with no finite upper cutoff.
   Physical validated range is unknown, not infinite. Input q² is GeV² and
   b is GeV^-1. At equal y, changing b preserves W; at fixed q, it changes W.
6. Datasets are independent of isotope mass and target fractions, immutable,
   and package files are checked against pinned hashes before loading.

## Checks performed

- Existing baseline: **730 passed** with warnings as errors.
- New focused regression tests: **100 passed**.
- Full suite: **830 passed** using `python -m pytest -q -W error`.
- Two independent exports with pinned SymPy produce identical bytes for every
  JSON file and the manifest.
- Original-source Mathics replay: 540 distinct operator/orbit/rank/coordinate
  values; Python transcription maximum absolute discrepancy about 1.34e-14.
  All 320 assembled W references agree; maximum scaled relative discrepancy
  about 4.2e-13. Tests use explicit machine-arithmetic tolerances.
- 465 printed coefficient checks from Fitzpatrick Appendix A.3 agree with the
  corresponding rounding intervals, including asymmetric mixed pairs. Printed
  polynomials are not high-precision finite-y references.
- Wheel packaging includes all six JSON files; the runtime loader has no
  SymPy/Mathics/Mathematica dependency.

The replay interpreter shares SymPy angular-symbol machinery with the exact
port. A tester with commercial Mathematica can add a valuable third check:
load the checksum-pinned original supplement, disable Helm, select each Ge
isotope, evaluate `ResponseNuclear[y,i,t,u]` at y=1/4 and 1, multiply by four,
and transpose t/u for i=7,8. Account explicitly for the source's `MyChop`
behavior: this export keeps tiny analytic terms instead of chopping them.
This extra interpreter check has not been claimed as completed here.

## Reproduction and evidence

Follow the pinned development commands in `GE_NUCLEAR_DATA.md`. Reference
fixtures are `tests/fixtures/ge_afh2014_source_replay.json` and
`ge_fitzpatrick_printed.json`; the latter retains source p/n mixed order.
The original supplement must pass its full-file SHA-256 gate before either
regeneration tool runs. Regenerated analytic data should match byte for byte;
source replay values should match within test tolerances across platforms.

Source SHA-256:
`943a403e216633840887bbfcbec55c56e0c5813b86335b49868acb89e74a49b9`.
Manifest SHA-256:
`b17df65326c24a46a4effef6f2b5881304c9c0f03792417c2635e921a6a76f33`.

| Artifact | SHA-256 |
|---|---|
| Ge-70 | `01c5a505d3dc6b2b0fe93e4a0195f2ff6dc95553af462f32d24211ac80ee12a5` |
| Ge-72 | `baba554841bd646a3034729214a3bee7c26f7a466262efc4b2a4f6362c8d6379` |
| Ge-73 | `7c97e13c266a751961bed8214512013e2d5941d1fbc8ada252de5248dc354256` |
| Ge-74 | `ddf30673b689a41ccce308c15d0e70a9fcf36820c68fed3e84ceb2eb232dee64` |
| Ge-76 | `d803be8ef578c06abb4519c36d1085fad9fabc4ab91bafe5bffbfe6a3169be20` |

The manifest records both exact generator script hashes. Each isotope JSON
also records the primary-source and supporting readme/audit hashes. The final
delivery includes full test output, exact commit/tree identity and a unified
diff against the baseline.

## Remaining source/model questions

- The particle-number closure deficit's origin is unresolved; no repair is
  authorized or applied.
- No empirical validity interval is established. The source replay and printed
  checks validate transcription, not the physical accuracy of the nuclear model
  at arbitrary q.
- Phi-tilde' is checked against source algebra/replay, not a separate printed
  Fitzpatrick table or alternative nuclear calculation.
- The GCN2850 label is supported by the audit and numerical fingerprint, not an
  explicit label in the primary file; original shell-model run artifacts and
  original APS download date were not provided.
- Analyses must explicitly select b, nuclear masses and isotope fractions. The
  supplement prescription is available as an explicit helper, not a hidden
  generic target default. These choices need their own numerical regressions.

No detector/statistical module, historical notebook, frozen convention or
production contraction was changed. No physical decay background, alternative
Ge nuclear dataset or coefficient renormalization was introduced.
