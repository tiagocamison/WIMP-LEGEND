# RATE-AUDIT-001 resolution: Phase A of GENERIC-ANALYSIS-001

Starting source: `44c4f2523b83b926a8b11d153ebfd5d8d04531a5`, exact
`origin/refactor` tip on 2026-10-09; parent
`efd56745fa5e35461122ae36afcfe105115061b1`. Work is isolated on
`generic-analysis-001`. No main/refactor update or push is performed.

Before edits: 492 repository tests pass with `-W error`. Recovered independent
scripts give 154 old-audit passes and 236 rate-audit passes / 1 failure. The
failure is exactly RATE-AUDIT-001: expected
`2.669722691466525e-291 events/(kg day keV)`, actual zero.

The final product now decomposes each finite factor with `math.frexp`,
accumulates a bounded mantissa and integer binary exponent, includes division
by the WIMP mass in the same representation, and restores scale once with
`math.ldexp`. No physical factor is changed: speed-flux integral, (hbar c)^2,
1e-6, c in cm/s, density/mass, nuclei/kg, and seconds/day. Exact zeros stay zero;
signed contributions stay signed. Overflow raises ValueError. Nonzero results
rounding below the smallest subnormal to zero also raise ValueError;
representable subnormals retain the usual reduced relative precision. This is
floating-point arithmetic, not arbitrary precision, and cannot recover prior
underflow in R/W evaluation, callback products, quadrature, or cancellations.

After edits: **504 repository tests pass**; **154 old independent audit tests
pass unchanged**; **237 new rate-audit tests pass unchanged**. The original
reproducer agrees with its independent 90-digit Decimal expectation (relative
tolerance 1e-12; zero absolute tolerance). All three recovered script SHA-256
hashes match the audit reports:

- test_independent.py: 1f1827d261b49a0928db7eec5c33c0bc5d518d976d107f245fd4de6639ac23e6
- test_falsifications.py: 99438a0db818631638f99d80b7282c8f6ed8d56505443abf4057b4e0005de6cb
- test_rate_audit.py: a90623aca9133d9293181cad00c4db5dbbb5fbc0f72cc72793612c172ba24aae

The new permanent `tests/test_rate_product.py` contains the verbatim original
reproducer and its synthetic fixture, ordinary/tiny-rate scaling and sign
checks, Decimal product references spanning intermediate overflow/underflow,
zero/subnormal cases, and genuine range-error checks: 12 added cases. Existing
tests and tolerances are unchanged. Runtime Python 3.12.14, NumPy 2.3.5,
SciPy 1.17.0, pytest 9.1.1. pytest was installed in a temporary dependency path.

Commands (from checkout, then workspace):

```bash
PYTHONPATH=/tmp/wimp-legend-test-deps python -m pytest -q -W error
PYTHONPATH=/tmp/wimp-legend-test-deps:WIMP-LEGEND/src python -m pytest -q -W error audit-work/test_independent.py audit-work/test_falsifications.py audit-work/test_rate_audit.py
```

This record is coder verification, not a new independent audit. The green
Phase-A tree is committed separately before any Phase-B implementation.
