"""FASTQ input validation (REQ-001, REQ-002).

The pipeline refuses to start on inputs that would otherwise fail late,
or worse, succeed silently with wrong results.
"""

from __future__ import annotations

import gzip
from dataclasses import dataclass
from pathlib import Path


class InputValidationError(ValueError):
    """Raised when an input file is unusable. The message names the file."""


@dataclass
class FastqSummary:
    path: str
    n_reads: int
    min_len: int
    max_len: int


def _open(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path, "rt")


def validate_fastq(path: str | Path, max_reads: int | None = None) -> FastqSummary:
    """Check a FASTQ file is well formed and return a short summary.

    Rejects: missing file, empty file, corrupt gzip, records that are not 4
    lines, header not starting with '@', separator not starting with '+',
    sequence/quality length mismatch, illegal bases.
    """
    path = Path(path)
    if not path.exists():
        raise InputValidationError(f"{path}: file not found")
    if path.stat().st_size == 0:
        raise InputValidationError(f"{path}: file is empty")

    allowed = set("ACGTN")
    n = 0
    min_len, max_len = 10**9, 0
    try:
        with _open(path) as fh:
            while True:
                header = fh.readline()
                if not header:
                    break
                seq = fh.readline().rstrip("\n")
                plus = fh.readline()
                qual = fh.readline().rstrip("\n")
                rec = n + 1
                if not seq or not plus or not qual:
                    raise InputValidationError(
                        f"{path}: truncated record {rec} (file ends mid-record)"
                    )
                if not header.startswith("@"):
                    raise InputValidationError(f"{path}: record {rec} header does not start with '@'")
                if not plus.startswith("+"):
                    raise InputValidationError(f"{path}: record {rec} separator line does not start with '+'")
                if len(seq) != len(qual):
                    raise InputValidationError(
                        f"{path}: record {rec} sequence length {len(seq)} != quality length {len(qual)}"
                    )
                bad = set(seq.upper()) - allowed
                if bad:
                    raise InputValidationError(f"{path}: record {rec} has illegal bases {sorted(bad)}")
                n += 1
                min_len = min(min_len, len(seq))
                max_len = max(max_len, len(seq))
                if max_reads and n >= max_reads:
                    break
    except (OSError, EOFError) as exc:  # corrupt or truncated gzip
        raise InputValidationError(f"{path}: cannot be read ({exc})") from exc

    if n == 0:
        raise InputValidationError(f"{path}: contains no reads")
    return FastqSummary(str(path), n, min_len, max_len)


def validate_pair(r1: str | Path, r2: str | Path) -> tuple[FastqSummary, FastqSummary]:
    """Validate both mates and check they have the same number of reads."""
    s1, s2 = validate_fastq(r1), validate_fastq(r2)
    if s1.n_reads != s2.n_reads:
        raise InputValidationError(
            f"{r1} has {s1.n_reads} reads but {r2} has {s2.n_reads}; mates are not paired"
        )
    return s1, s2
