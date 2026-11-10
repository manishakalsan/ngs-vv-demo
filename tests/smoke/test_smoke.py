"""Smoke tests: the build is alive enough to be worth testing further (REQ-001, REQ-010).

Run first in CI. If these fail, nothing else is meaningful.
"""

import gzip

import pytest


@pytest.mark.req("REQ-001")
def test_pipeline_completes_and_produces_indexed_vcf(baseline_run):
    out, s = baseline_run["out"], baseline_run["sample"]
    vcf = out / s / f"{s}.vcf.gz"
    assert vcf.exists() and vcf.stat().st_size > 0
    assert (out / s / f"{s}.vcf.gz.tbi").exists()
    assert (out / s / f"{s}.bam.bai").exists()
    with gzip.open(vcf, "rt") as fh:
        records = [l for l in fh if not l.startswith("#")]
    assert len(records) > 0, "VCF has no variant records"


@pytest.mark.req("REQ-001")
def test_qc_reports_written(baseline_run):
    qc = baseline_run["out"] / "SYNTH1" / "qc"
    for name in ("input_validation.json", "SYNTH1.fastp.json", "SYNTH1.fastp.html"):
        assert (qc / name).exists(), name


@pytest.mark.req("REQ-010")
def test_software_versions_recorded(baseline_run):
    text = (baseline_run["out"] / "software_versions.txt").read_text()
    for tool in ("ngs-vv-demo", "snakemake", "fastp", "bwa", "samtools", "bcftools"):
        assert tool in text, f"{tool} version not recorded"
    # every line carries a version number
    assert all(any(ch.isdigit() for ch in line) for line in text.strip().splitlines())
