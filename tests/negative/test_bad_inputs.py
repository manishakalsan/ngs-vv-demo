"""Negative tests at pipeline level (REQ-002, REQ-003).

Each bad input must: exit non-zero, print an error naming the problem,
and stop BEFORE alignment (no BAM written). Failing late, or worse
succeeding, is the defect these tests guard against.
"""

import gzip
import shutil

import pytest

from conftest import run_pipeline


def assert_failed_before_alignment(proc, workdir, expect_in_stderr):
    assert proc.returncode != 0, "pipeline should have failed but succeeded"
    assert expect_in_stderr in proc.stderr, proc.stderr[-2000:]
    assert not list((workdir / "results").rglob("*.bam")), "BAM was written: did not fail fast"


def fastq_override(r1, r2):
    return {"samples": {"SYNTH1": {"r1": str(r1), "r2": str(r2)}}}


@pytest.mark.req("REQ-002")
def test_truncated_fastq_fails_fast(tmp_path, synthetic, tools_available):
    bad = tmp_path / "trunc_R1.fastq.gz"
    with gzip.open(synthetic["r1"], "rt") as fh:
        lines = fh.readlines()
    with gzip.open(bad, "wt") as fh:
        fh.writelines(lines[:-2])  # last record loses its '+' and quality lines
    proc = run_pipeline(tmp_path / "wd", fastq_override(bad, synthetic["r2"]))
    assert_failed_before_alignment(proc, tmp_path / "wd", "trunc_R1.fastq.gz")


@pytest.mark.req("REQ-002")
def test_corrupt_gzip_fails_fast(tmp_path, synthetic, tools_available):
    bad = tmp_path / "corrupt_R1.fastq.gz"
    data = synthetic["r1"].read_bytes()
    bad.write_bytes(data[: len(data) // 3])  # interrupted download
    proc = run_pipeline(tmp_path / "wd", fastq_override(bad, synthetic["r2"]))
    assert_failed_before_alignment(proc, tmp_path / "wd", "corrupt_R1.fastq.gz")


@pytest.mark.req("REQ-002")
def test_empty_fastq_fails_fast(tmp_path, synthetic, tools_available):
    bad = tmp_path / "empty_R1.fastq.gz"
    bad.write_bytes(b"")
    proc = run_pipeline(tmp_path / "wd", fastq_override(bad, synthetic["r2"]))
    assert_failed_before_alignment(proc, tmp_path / "wd", "empty")


@pytest.mark.req("REQ-002")
def test_mates_from_different_samples_fail_fast(tmp_path, synthetic, tools_available):
    r2 = tmp_path / "short_R2.fastq.gz"
    with gzip.open(synthetic["r2"], "rt") as fh:
        lines = fh.readlines()
    with gzip.open(r2, "wt") as fh:
        fh.writelines(lines[:4000])  # 1,000 reads instead of 6,000
    proc = run_pipeline(tmp_path / "wd", fastq_override(synthetic["r1"], r2))
    assert_failed_before_alignment(proc, tmp_path / "wd", "not paired")


@pytest.mark.req("REQ-002")
def test_missing_fastq_fails(tmp_path, synthetic, tools_available):
    proc = run_pipeline(tmp_path / "wd", fastq_override(tmp_path / "nope_R1.fastq.gz", synthetic["r2"]))
    assert proc.returncode != 0
    assert "nope_R1.fastq.gz" in proc.stderr + proc.stdout


@pytest.mark.req("REQ-003")
def test_wrong_contig_naming_fails_fast(tmp_path, synthetic, tools_available):
    """Reference uses 'chr' naming the run does not expect (e.g. hg19 vs GRCh37)."""
    ref = tmp_path / "ref_renamed.fa"
    ref.write_text(synthetic["ref"].read_text().replace(">synth1", ">chrsynth1"))
    proc = run_pipeline(tmp_path / "wd", {"reference": str(ref)})
    assert_failed_before_alignment(proc, tmp_path / "wd", "expected contig 'synth1' is missing")


@pytest.mark.req("REQ-003")
def test_wrong_build_length_fails_fast(tmp_path, synthetic, tools_available):
    """Same contig name, different sequence length = a different assembly."""
    ref = tmp_path / "ref_short.fa"
    lines = synthetic["ref"].read_text().splitlines()
    ref.write_text("\n".join(lines[:-10]) + "\n")
    proc = run_pipeline(tmp_path / "wd", {"reference": str(ref)})
    assert_failed_before_alignment(proc, tmp_path / "wd", "wrong genome build")


@pytest.mark.req("REQ-003")
def test_bed_on_other_build_fails_fast(tmp_path, synthetic, tools_available):
    bed = tmp_path / "confident_chr.bed"
    bed.write_text("chr20\t500\t29500\n")
    proc = run_pipeline(tmp_path / "wd", {"confident_bed": str(bed)})
    assert_failed_before_alignment(proc, tmp_path / "wd", "not in reference")
