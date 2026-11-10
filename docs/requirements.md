# Software requirements specification

**Product:** ngs-vv-demo germline small-variant pipeline · **Version:** 1.0.0 · **Status:** demonstration, not a regulated product

The machine-readable list is [`requirements.yaml`](requirements.yaml); it is the single source the traceability matrix is generated from. This page explains each requirement and why it exists.

## Intended use (demo)

Detect germline SNVs and small indels (≤ 50 bp) from paired-end Illumina short reads against a single reference contig or region. Research-use demonstration only; no clinical claims.

## Requirements

| ID | Requirement | Why it exists | Verified by |
| --- | --- | --- | --- |
| REQ-001 | Accept paired-end gzipped FASTQ; produce an indexed VCF | Core function | smoke, unit |
| REQ-002 | Reject malformed, truncated, empty, missing or unpaired FASTQ before alignment; error names the file | An interrupted download or a swapped mate file can otherwise produce a plausible but wrong VCF | unit, negative |
| REQ-003 | Reject a reference whose contigs do not match the expected build or region BED, before alignment | Build mismatch (GRCh37 vs GRCh38, `20` vs `chr20`) is a classic silent failure | unit, negative |
| REQ-004 | SNV PPA ≥ 0.99 and PPV ≥ 0.99 in confident regions | Performance claim | regression |
| REQ-005 | Indel PPA ≥ 0.95 and PPV ≥ 0.95 in confident regions | Indels are harder; a separate claim stops a good SNV score hiding poor indel calling | regression |
| REQ-006 | Identical input and config give identical variant records across repeated runs and thread counts | Repeatability; multithreaded aligners can be non-deterministic | regression |
| REQ-007 | Test dataset completes in ≤ 120 s on 2 cores | Performance baseline for trend tracking | performance |
| REQ-008 | Every critical variant is called PASS with the correct genotype | Stand-in for "clinically actionable variants must never be missed" | regression |
| REQ-009 | Output VCF sorted, normalised (left-aligned, biallelic); QUAL < threshold marked LowQual | Downstream tools and truth comparison depend on canonical representation | unit, regression |
| REQ-010 | Each run records versions of all software components | Configuration management; reproduce any result | smoke |
| REQ-011 | Evaluation code computes TP/FP/FN/PPA/PPV correctly | The test tool itself must be validated, or every accuracy claim is unsupported | unit, regression |
| REQ-012 | Variants outside confident regions excluded from metrics | Truth is only reliable inside confident regions | unit |

## How acceptance thresholds were set

Thresholds come from the intended use, then are confirmed on the v1.0.0 baseline run and **locked** in `tests/regression/acceptance_criteria.yaml`. They are not adjusted to match whatever a later release produces. Changing them requires a change request (see [test plan](test_plan.md#6-change-control)).

On the synthetic dataset v1.0.0 scores 1.0 on all metrics, because the data are easy by design. The thresholds are set for realistic data, which is why the [GIAB benchmark](giab_mode.md) matters.
