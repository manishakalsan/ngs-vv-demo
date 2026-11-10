# SOUP list (software of unknown provenance)

IEC 62304 calls third-party software not developed under the manufacturer's own process "SOUP". Each item needs its version, the function it performs, and known anomalies that could affect safety. All versions are pinned in `envs/environment.yml` and, byte-exact, in `envs/conda-linux-64.lock`; each run also writes `software_versions.txt`.

| Component | Version | Function in pipeline | Known issues relevant here | Risk control |
| --- | --- | --- | --- | --- |
| fastp | 0.23.4 | Adapter/quality trimming | Multithreaded output order can vary | Run with `-w 1`; repeatability test (REQ-006) |
| BWA-MEM | 0.7.18 | Read alignment | Output depends on thread count unless batch size is fixed | `-K 10000000`; REQ-006 test runs 1, 2 and 4 threads |
| samtools | 1.21 | Sort, index, FASTA index | None identified for this use | Smoke tests |
| bcftools mpileup/call | 1.21 | Variant calling | Weaker indel calling than haplotype-based callers; `-I` silently skips indels | Separate indel threshold (REQ-005); CR-001 shows `-I` is caught |
| bcftools norm | 1.21 | Left-align and split multi-allelics | Requires matching reference | REQ-003 build check; REQ-009 re-normalisation test |
| bcftools filter | 1.21 | Mark QUAL < 20 as LowQual | Expression errors fail silently to no-op if mistyped | REQ-009 test checks every record's FILTER against its QUAL |
| Snakemake | 8.30.0 | Workflow engine | Config files merge recursively (dicts merge, lists replace) | `expected_contigs` is a list; `run_samples` selector |
| Python | 3.11 | Validation and evaluation code | — | Unit tests, 99% coverage |
| PyYAML | 6.0.3 | Config parsing | — | — |
| pytest / pytest-cov | 9.2.0 / 7.1.0 | Test execution (not shipped) | — | Tool used for verification only |

Anomaly lists for each tool should be reviewed against its release notes/issue tracker at every version change; this demo has not done that formal review.
