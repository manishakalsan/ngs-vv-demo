"""Unit tests: genome-build guard (REQ-003)."""

import pytest

from vvdemo.fastq import InputValidationError
from vvdemo.reference import bed_contigs, check_reference, read_fai, read_fasta_contigs


@pytest.fixture
def ref(tmp_path):
    p = tmp_path / "ref.fa"
    p.write_text(">chr20 some description\nACGTACGTAC\nACGT\n>chr21\nAAAA\n")
    return p


@pytest.mark.req("REQ-003")
def test_contig_lengths_read_from_fasta(ref):
    assert read_fasta_contigs(ref) == {"chr20": 14, "chr21": 4}


@pytest.mark.req("REQ-003")
def test_matching_reference_passes(ref):
    assert "chr20" in check_reference(ref, {"chr20": 14})


@pytest.mark.req("REQ-003")
@pytest.mark.parametrize("expected, message", [
    ({"20": 14}, "missing"),          # naming mismatch: 20 vs chr20 (GRCh37 vs hg19 style)
    ({"chr20": 99}, "length"),        # same name, different build
])
def test_wrong_build_rejected(ref, expected, message):
    with pytest.raises(InputValidationError, match=message):
        check_reference(ref, expected)


@pytest.mark.req("REQ-003")
def test_bed_on_other_contig_naming_rejected(ref, tmp_path):
    bed = tmp_path / "regions.bed"
    bed.write_text("track name=x\n20\t0\t10\n")
    with pytest.raises(InputValidationError, match="not in reference"):
        check_reference(ref, None, bed)
    assert bed_contigs(bed) == {"20"}


@pytest.mark.req("REQ-003")
def test_missing_or_empty_reference_rejected(tmp_path):
    with pytest.raises(InputValidationError, match="not found"):
        check_reference(tmp_path / "none.fa")
    empty = tmp_path / "empty.fa"
    empty.write_text("")
    with pytest.raises(InputValidationError, match="no FASTA records"):
        check_reference(empty)


@pytest.mark.req("REQ-003")
def test_read_fai(tmp_path):
    fai = tmp_path / "ref.fa.fai"
    fai.write_text("chr20\t14\t7\t10\t11\n\n")
    assert read_fai(fai) == {"chr20": 14}
