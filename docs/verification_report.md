# Verification report: ngs-vv-demo 1.0.0

**Report ID:** VR-001 · **Date:** 2026-11-10 · **Plan:** [TP-001](test_plan.md) · **Environment:** `envs/conda-linux-64.lock` (linux-64), 2 cores

## Conclusion

All 71 tests passed and all 12 requirements are covered and passing; the release gate passed. One open risk (H-09, difficult genomic regions) is outside what the synthetic CI data can verify and is carried forward.

## Results by level

| Level | Tests | Passed |
| --- | --: | --: |
| Unit | 47 | 47 |
| Smoke | 3 | 3 |
| Negative | 8 | 8 |
| Regression | 12 | 12 |
| Performance | 1 | 1 |
| **Total** | **71** | **71** |

Unit coverage of `src/vvdemo`: 100% of lines, 98.8% including branches (gate: 90%).

## Accuracy on the synthetic dataset (confident regions)

| Class | Truth | TP | FP | FN | PPA | PPV | Threshold met |
| --- | --: | --: | --: | --: | --: | --: | --- |
| SNV | 40 | 40 | 0 | 0 | 1.000 | 1.000 | Yes (≥ 0.99 / ≥ 0.99) |
| Indel | 12 | 12 | 0 | 0 | 1.000 | 1.000 | Yes (≥ 0.95 / ≥ 0.95) |

All 4 critical variants called PASS with correct genotype. The SNV at position 200, outside the confident BED, was called and correctly excluded from metrics. These scores reflect easy synthetic data and are **not** a performance claim for real samples.

## Other results

| Check | Result |
| --- | --- |
| Repeatability, 1 / 2 / 4 threads (REQ-006) | VCF records identical (MD5 of body) |
| Runtime (REQ-007) | 1.8 s wall vs 120 s limit |
| Output format (REQ-009) | Sorted, biallelic; `bcftools norm -c e` re-run changes 0 records; 53/53 records PASS with QUAL ≥ 20 |
| Versions recorded (REQ-010) | ngs-vv-demo 1.0.0, Snakemake 8.30.0, fastp 0.23.4, BWA 0.7.18, samtools 1.21, bcftools 1.21 |

## Anomalies and deviations found during verification

| # | Found by | Issue | Resolution |
| --- | --- | --- | --- |
| 1 | Coverage report | A gzipped empty FASTQ (~20 bytes) passed the "file is empty" size check | Already rejected by the "contains no reads" branch, but untested; test added |
| 2 | Regression test | `bcftools norm` summary format differs between versions (a `removed` field was added) | Parser made field-name based |
| 3 | Pipeline review | Snakemake merges `--configfile` dicts, so the GIAB config would also demand the synthetic contig | `expected_contigs` changed to a list; `run_samples` added |

## Change impact: CR-001 (rejected)

Proposed: skip indel calling to speed up `mpileup`. The regression suite failed 4 tests: indel PPA fell from 1.0 to 0.0, and the critical-variant check failed. Full record: [changes/CR-001_impact.md](changes/CR-001_impact.md).

## Not verified here

- GIAB benchmark: scripts and config provided; not run in this report because the data could not be downloaded in the build environment.
- Docker image: Dockerfile provided; the lock file it installs from was verified to recreate the environment and pass the tests, but the image itself was not built.
