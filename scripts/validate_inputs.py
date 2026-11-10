#!/usr/bin/env python3
"""Pre-flight input validation. First rule of the pipeline.

Exits non-zero with a message naming the bad file, so the run stops before
alignment (REQ-002, REQ-003). Writes a JSON summary on success.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from vvdemo.fastq import InputValidationError, validate_pair  # noqa: E402
from vvdemo.reference import check_reference  # noqa: E402


def parse_expected(items: list[str]) -> dict[str, int]:
    out = {}
    for item in items or []:
        name, _, length = item.partition(":")
        out[name] = int(length) if length else 0
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--r1", required=True)
    ap.add_argument("--r2", required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--bed", action="append", default=[], help="region BEDs that must match the reference")
    ap.add_argument("--expected-contig", action="append", default=[], help="NAME[:LENGTH]")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    try:
        s1, s2 = validate_pair(a.r1, a.r2)
        contigs = check_reference(a.reference, parse_expected(a.expected_contig))
        for bed in a.bed:
            check_reference(a.reference, None, bed)
    except InputValidationError as exc:
        print(f"INPUT VALIDATION FAILED: {exc}", file=sys.stderr)
        return 2

    summary = {
        "status": "PASS",
        "r1": vars(s1), "r2": vars(s2),
        "reference": a.reference,
        "contigs": contigs,
    }
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text(json.dumps(summary, indent=2))
    print(f"input validation PASS: {s1.n_reads} read pairs, contigs {list(contigs)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
