# Risk analysis (ISO 14971-style, simplified)

Hazard = a way the software could contribute to a wrong result. Severity and likelihood on a 1–3 scale before controls; residual risk after controls are verified.

| ID | Hazard | Cause | Sev | Lik | Control | Verified by | Residual |
| --- | --- | --- | :-: | :-: | --- | --- | --- |
| H-01 | Variants missed or invented without warning | Truncated/corrupt FASTQ processed as if complete | 3 | 2 | Pre-flight FASTQ validation, fail before alignment | REQ-002 tests | Low |
| H-02 | Mates from different samples combined | Sample-sheet error | 3 | 2 | Mate read-count check | `test_mates_from_different_samples_fail_fast` | Low (count check cannot catch equal-sized swaps; read-name check would) |
| H-03 | Coordinates wrong for downstream annotation | Wrong genome build or contig naming | 3 | 2 | Contig name + length check against expected build | REQ-003 tests | Low |
| H-04 | Indels not reported | Tool option or version change | 3 | 2 | Separate indel threshold; critical-variant list; regression job | REQ-005, REQ-008, CR-001 | Low |
| H-05 | Different results on re-analysis | Non-deterministic multithreading | 2 | 2 | Fixed BWA batch size, single-thread fastp | REQ-006 | Low |
| H-06 | Accuracy claims overstated | Bug in evaluation code | 3 | 1 | Hand-counted unit tests of the comparator; truth normalisation check | REQ-011 | Low |
| H-07 | Low-confidence calls reported as confident | Filter not applied | 2 | 2 | QUAL filter + per-record FILTER check | REQ-009 | Low |
| H-08 | Results cannot be reproduced later | Unrecorded tool versions | 2 | 2 | Pinned env, lock file, versions written per run | REQ-010 | Low |
| H-09 | Poor accuracy in difficult regions | Repeats, homopolymers, low mappability | 3 | 3 | Not controlled by synthetic data | GIAB benchmark (manual) | **Open**: needs stratified GIAB analysis |

H-09 is the honest gap: the CI dataset cannot show it. That is the next piece of work.
