# CODE_AUDIT_001 resolution — Phase A

Audited base: `efd56745fa5e35461122ae36afcfe105115061b1`, fetched from
`origin/refactor`. Initial local HEAD was `5d5b22724f746b75a73e0931aeffa0492f2abe1f`;
its complete tree matches the audited commit. Work is isolated on
`audit-remediation-reference-rates`, based on the exact audited SHA.
The independent report is unchanged. This is a coder remediation record, not a
new independent certification.

## Reproduction and gate

Before production changes: **346 repository tests passed**; the independent
Appendix A scripts yielded **144 passed, 10 failed**, exactly the reported cases.
Extracted script SHA-256 hashes match Appendix C:

- `test_independent.py`: `1f1827d261b49a0928db7eec5c33c0bc5d518d976d107f245fd4de6639ac23e6`
- `test_falsifications.py`: `99438a0db818631638f99d80b7282c8f6ed8d56505443abf4057b4e0005de6cb`

After remediation: **383 repository tests passed**, including 37 focused new
cases; **154/154 independent audit tests passed unchanged**. No skips or warnings.
No pre-existing tests or tolerances were changed. Both commands use `-q -W error`.
The external scripts remain outside the checkout in `../audit-work`, with the
checkout named `audit-repo` as required by their import-boundary test.

```bash
PYTHONPATH=/tmp/wimp-legend-test-deps python -m pytest -q -W error
# From ../audit-work:
PYTHONPATH=/tmp/wimp-legend-test-deps:../audit-repo/src python -m pytest test_independent.py test_falsifications.py -q -W error
```

Runtime: Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0, pytest 9.1.1. Temporary
pytest installation is environmental; project dependencies are unchanged.

## Findings and permanent regressions

All test names below are in `tests/test_audit_regressions.py`.

| Finding | Production change | Permanent regression | Original reproducer |
| --- | --- | --- | --- |
| AUDIT-001 | Validate kind, transformation, fraction_kind and momentum variable as scalar strings before equality/membership; reject arrays. Store accepted scalar text as Python strings, preserving exact content. Also review/copy reference, locator, notes, citations, symbol, channel, metadata fields, labels and dataset ID. | `test_audit001_source_text_rejects_mutable_arrays`, `test_audit001_metadata_text_rejects_mutable_arrays`, `test_audit001_momentum_and_fraction_tags`, `test_audit001_other_text_and_exact_owned_ids` | All four pass |
| AUDIT-002 | Factor only the Sigma-double-prime c4/c6 bilinear as `(a4+Q*a6)*(b4+Q*b6)`. No physics change or clipping. | `test_audit002_factored_sigma_high_precision`: 80-digit Decimal, both sides of cancellation, exact zero, all ordered pairs, negative off-diagonal values | Pass |
| AUDIT-003 | Compute q from square roots before multiplication, preserving the exact 1e-6 conversion. Use mantissa/exponent arithmetic for q/(2 mu); reject nonfinite or erased nonzero outputs. | `test_audit003_scaled_momentum_and_speed`, `test_audit003_nonrepresentable_speed_rejected`: Decimal, extreme and ordinary scales | Both pass |
| AUDIT-004 | Evaluate the full/partial/zero-boost PDFs as a dimensionless shape divided once by v0; check finite density output. | `test_audit004_dimensionless_pdf_scaling_and_normalization`: ordinary-to-tiny scale covariance, normalization, zero/tiny/equal/above-escape boosts | Pass |
| AUDIT-005 | Table coordinates are x=v/vmax; weight the scaled density vmax*f(vmax*x). Integrate the first cell adaptively in x, then apply the dimensional moment scale. Check weighted samples, cumulative values and final table; reject nonrepresentable tables. Clamp out-of-support thresholds before scaling. | `test_audit005_dimensionless_table`, `test_audit005_unrepresentable_table_fails`: analytic polynomial moments, breakpoints, shapes, explicit range failure | Pass |
| AUDIT-006 | Range-aware average: sum before halving when both operands are at most half maximum float magnitude; otherwise halve before adding. Reuse for difference with negated second input. | `test_audit006_equal_and_opposite_range`: smallest subnormal through 1e308, same and opposite signs | Pass |

The existing suites retain ordinary SHM normalization/support, moments and table
convergence, kinematic inversion, composition conversions/provenance, enforced
source ranges, large Wilson conversions, ordered interference, randomized
Eq. (38) checks, and import boundaries. The independent audit adds its original
Decimal, angular, extended-precision and spin-operator references.

## Retained limits

These fixes do not claim arbitrary-precision support over every possible finite
input combination. The SHM constructor still requires a representable positive
3D normalization. Extreme speed ratios can exceed intermediate floating-point
range and fail explicitly; nonfinite density outputs are rejected. Moment tables
retain power >= -1 and O(v²) PDF assumptions, finite-dimensional-scale limits,
and linear-interpolation tail error. Unsupported moment scales fail explicitly;
reference adaptive near-singular quadrature can still raise IntegrationWarning.
True sub-float results need not retain arbitrary relative precision. Very large
Wilson inverse conversions and responses can remain nonrepresentable and fail.
No clipping or tolerance relaxation was used.

AUDIT-007 remains a historical-document note, and AUDIT-008 remains an independent
audit-procedure limitation; neither warrants a production physics change.
The frozen normative document and convention ID are unchanged.

Phase B may begin only after this green checkpoint is committed.
