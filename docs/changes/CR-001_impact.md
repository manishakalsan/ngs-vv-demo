# CR-001 impact assessment: skip indel calling

Change: `caller.call_indels: false` (proposed to speed up `bcftools mpileup`).

| Class | Metric | v1.0.0 (released) | v1.1.0-rc1 (CR-001) | Change |
| --- | --- | --: | --: | --: |
| SNV | TP | 40 | 40 | +0 |
| SNV | FN | 0 | 0 | +0 |
| SNV | FP | 0 | 0 | +0 |
| SNV | PPA | 1.0 | 1.0 | +0 |
| SNV | PPV | 1.0 | 1.0 | +0 |
| INDEL | TP | 12 | 0 | -12 |
| INDEL | FN | 0 | 12 | +12 |
| INDEL | FP | 0 | 0 | +0 |
| INDEL | PPA | 1.0 | 0.0 | -1 |
| INDEL | PPV | 1.0 | None |  |

Regression tests failed: 4

- `tests/regression/test_accuracy.py::test_indel_accuracy_meets_threshold`
- `tests/regression/test_accuracy.py::test_no_drift_from_locked_baseline[PPA-INDEL]`
- `tests/regression/test_accuracy.py::test_no_drift_from_locked_baseline[PPV-INDEL]`
- `tests/regression/test_accuracy.py::test_critical_variants_always_called`

**Decision: REJECT.** REQ-005 (indel accuracy) and REQ-008 (critical variants) are violated.
