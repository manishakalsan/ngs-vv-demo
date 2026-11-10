"""Unit tests: FASTQ validation (REQ-002). Positive and negative cases."""

import gzip

import pytest

from vvdemo.fastq import InputValidationError, validate_fastq, validate_pair

GOOD = "@r1\nACGTN\n+\nFFFFF\n@r2\nACG\n+\nFFF\n"


def write(tmp_path, name, text, gz=False):
    p = tmp_path / name
    if gz:
        with gzip.open(p, "wt") as fh:
            fh.write(text)
    else:
        p.write_text(text)
    return p


@pytest.mark.req("REQ-001")
@pytest.mark.parametrize("gz", [False, True])
def test_valid_fastq_is_summarised(tmp_path, gz):
    p = write(tmp_path, "ok.fq.gz" if gz else "ok.fq", GOOD, gz)
    s = validate_fastq(p)
    assert (s.n_reads, s.min_len, s.max_len) == (2, 3, 5)


@pytest.mark.req("REQ-002")
@pytest.mark.parametrize("text, message", [
    ("@r1\nACGT\n+\nFFF\n", "quality length"),           # seq/qual length mismatch
    ("r1\nACGT\n+\nFFFF\n", "does not start with '@'"),   # bad header
    ("@r1\nACGT\n-\nFFFF\n", "separator"),                # bad separator
    ("@r1\nACXT\n+\nFFFF\n", "illegal bases"),            # non-ACGTN base
    ("@r1\nACGT\n+\nFFFF\n@r2\nACGT\n", "truncated"),     # file ends mid-record
])
def test_malformed_fastq_rejected(tmp_path, text, message):
    p = write(tmp_path, "bad.fq", text)
    with pytest.raises(InputValidationError, match=message) as exc:
        validate_fastq(p)
    assert "bad.fq" in str(exc.value), "error must name the file"


@pytest.mark.req("REQ-002")
def test_missing_file_rejected(tmp_path):
    with pytest.raises(InputValidationError, match="not found"):
        validate_fastq(tmp_path / "nope.fq.gz")


@pytest.mark.req("REQ-002")
def test_empty_file_rejected(tmp_path):
    p = write(tmp_path, "empty.fq", "")
    with pytest.raises(InputValidationError, match="empty"):
        validate_fastq(p)


@pytest.mark.req("REQ-002")
def test_corrupt_gzip_rejected(tmp_path):
    p = write(tmp_path, "good.fq.gz", GOOD * 50, gz=True)
    data = p.read_bytes()
    p.write_bytes(data[: len(data) // 2])  # truncated download
    with pytest.raises(InputValidationError, match="good.fq.gz"):
        validate_fastq(p)


@pytest.mark.req("REQ-002")
def test_unpaired_mates_rejected(tmp_path):
    r1 = write(tmp_path, "a_R1.fq", GOOD)
    r2 = write(tmp_path, "a_R2.fq", GOOD.split("@r2")[0])
    with pytest.raises(InputValidationError, match="not paired"):
        validate_pair(r1, r2)


@pytest.mark.req("REQ-001")
def test_max_reads_stops_early(tmp_path):
    p = write(tmp_path, "ok.fq", GOOD * 10)
    assert validate_fastq(p, max_reads=3).n_reads == 3


@pytest.mark.req("REQ-002")
def test_gzipped_empty_file_rejected(tmp_path):
    """Found via coverage: a .gz of nothing is ~20 bytes, so the size check alone misses it."""
    p = write(tmp_path, "nothing.fq.gz", "", gz=True)
    assert p.stat().st_size > 0
    with pytest.raises(InputValidationError, match="contains no reads"):
        validate_fastq(p)


@pytest.mark.req("REQ-001")
def test_valid_pair_accepted(tmp_path):
    s1, s2 = validate_pair(write(tmp_path, "a_R1.fq", GOOD), write(tmp_path, "a_R2.fq", GOOD))
    assert s1.n_reads == s2.n_reads == 2
