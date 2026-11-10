# Software verification test plan

**Product:** ngs-vv-demo 1.0.0 · **Plan ID:** TP-001

## 1. Scope

In scope: input validation, alignment-to-VCF pipeline, VCF output format, accuracy against truth, repeatability, runtime, evaluation tooling.
Out of scope (would be required for a real IVD): wet-lab steps, multiple instruments/reagent lots, clinical validation, UI, cybersecurity, user-documentation review.

## 2. Approach

Risk-based: tests concentrate where a failure would produce a **wrong result that looks right** (requirements marked risk H). Every test is linked to at least one requirement with `@pytest.mark.req(...)`; the traceability matrix is generated from those links, so it cannot drift from the code.

## 3. Test levels

| Level | Folder | What it proves | Runs |
| --- | --- | --- | --- |
| Unit | `tests/unit` | Each helper function is correct, incl. boundary and invalid values | every push |
| Smoke | `tests/smoke` | Pipeline runs end to end and writes all outputs | every push |
| Negative | `tests/negative` | Bad inputs fail fast, before alignment, with a clear error | every push |
| Regression | `tests/regression` | Accuracy, output format and repeatability unchanged vs locked baseline | every push and every release |
| Performance | `tests/performance` | Runtime within limit; recorded for trend | every push |
| Change impact | `scripts/demo_regression_catch.sh` | The suite detects a known-harmful change ("testing the tests") | every push (CI job 2) |

## 4. Test data

| Dataset | Purpose | Truth |
| --- | --- | --- |
| Synthetic, 30 kb, 6,000 pairs (2×100 bp, ~40×, 0.2% errors) | CI: fast, deterministic, exact truth | 41 SNVs + 12 indels inserted by `scripts/simulate_test_data.py`; 1 SNV deliberately outside the confident BED |
| GIAB HG002 GRCh38 chr20 slice | Realistic benchmark | GIAB v4.2.1 benchmark VCF + confident BED |
| Hand-written mini files in unit tests | Edge cases | Counted by hand |

## 5. Entry and exit criteria

Entry: environment created from `envs/` with pinned versions; unit tests pass.
Exit (release gate, enforced by `scripts/build_traceability.py`):
- every requirement has ≥ 1 linked test
- every linked test passes
- no test is unlinked
- unit coverage of `src/vvdemo` ≥ 90% (enforced by `pyproject.toml`)

## 6. Change control

Any change to code, config, tool versions or the reference is a change request (CR) with:
1. Description and reason.
2. Impact assessment: which outputs and requirements could be affected.
3. Run the full suite with the change applied (`VVDEMO_CONFIG_OVERLAY=<file> pytest`).
4. Compare against the locked baseline; investigate every changed call.
5. Decision (accept / reject) recorded in `docs/changes/`.

Worked example: [CR-001](changes/CR-001_impact.md), rejected.

## 7. Defect severity

| Severity | Definition | Example | Release? |
| --- | --- | --- | --- |
| Critical | Wrong or missing result without warning | Critical variant not called; bad FASTQ accepted | Blocks |
| Major | Function fails visibly, no workaround | Pipeline crashes on valid input | Blocks |
| Minor | Workaround exists | Unclear error text | Allowed if documented |
| Cosmetic | No functional effect | Typo in log | Allowed |

## 8. Known limitations of the test approach

- Truth comparison matches on normalised (CHROM, POS, REF, ALT); it is not haplotype-aware. Use hap.py for GIAB (see [giab_mode.md](giab_mode.md)).
- Synthetic data have no repeats, GC extremes or mapping ambiguity, so they cannot reveal problems in difficult regions.
- Runtime limit is loose because CI hardware varies.
