# Post-poster consolidation — 2026-11-18

`main` is the canonical development baseline again. It was fast-forwarded from
`72af6a1fc178bf97e241588ec404db8235c5a47e` to
`86cf0bd3beeb138cdc1e00db87992923720ace20`, preserving every commit and the user's
final poster plots. No physics code, datasets or frozen conventions were changed.
The baseline tree is `a9698ef83ed099c998df81822a41a042205f30f4`.

## Branch audit

Every original development tip below is an ancestor of `86cf0bd`. There were no
unmerged changes or conflicts between those tips. They do not need separate PRs.

| Original branch | Preserved commit | Disposition |
|---|---|---|
| main | 72af6a1fc178bf97e241588ec404db8235c5a47e | Historical archive |
| refactor | 93a8ed143578ad4b6a74bd636ada2943fa7e0ecb | Fully incorporated into main |
| poster | 874e59cd5cc45b923add3a3ebf10f26e1111e070 | Historical poster branch retained |
| poster-speedup | 43c063ecd24f8d8dac488481ca1d656f771f5fdb | Fully incorporated into main |
| poster-workbench | 406e3e3ded75fb92c1d2e58aa8ec6bb7a0f74678 | Fully incorporated into main |
| feature/pre-smearing-acceptance | 86cf0bd3beeb138cdc1e00db87992923720ace20 | Fully incorporated into main; final poster snapshot |

Permanent snapshot branches created before consolidation:

- `archive/2026-11-18/main`
- `archive/2026-11-18/poster`
- `archive/2026-11-18/poster-final`

Redundant original branch names remain on GitHub. The available authenticated
interface supports ref creation/update, but not deletion, and command-line Git
push is unauthenticated. No branch deletion is claimed. The commits remain
reachable through main and the snapshots. Once deletion access is available,
`refactor`, `poster-speedup`, `poster-workbench`, and
`feature/pre-smearing-acceptance` can be deleted after rechecking their tips have
not advanced. Keep the poster snapshot as a historical scientific reference.

## Verification

Current baseline: **278 passed in 33.83 s**, warnings treated as errors:

```sh
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python -m pytest -q -W error \
  tests/test_ionization_acceptance.py tests/test_ge_nuclear_data.py \
  tests/test_interactions.py tests/test_si_normalization.py
```

External pytest plugin auto-discovery was disabled because the environment's
unrelated AnyIO plugin raises a pytest deprecation warning during startup.
No project test was weakened or modified. This is focused current verification,
not a new full-suite run. Existing full-suite evidence (1140 passed) and the
later ionization checks remain historical evidence, not rerun claims.

Combined final focused run, including the five new RAPIDD artifact checks:
**283 passed in 32.52 s**, warnings treated as errors. See
[the comparison report](RAPIDD_COMPARISON.md) before using agreement with another
implementation as publication evidence.
