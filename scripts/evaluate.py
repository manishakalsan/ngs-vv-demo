#!/usr/bin/env python3
"""Compare pipeline calls to the truth set inside confident regions.

Writes metrics JSON (per variant class: TP, FP, FN, PPA, PPV) and a short
TSV. The regression tests read the JSON and apply the locked thresholds.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from vvdemo.vcf import compare, read_bed, read_vcf  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--calls", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--bed", required=True)
    ap.add_argument("--json", required=True)
    ap.add_argument("--tsv", required=True)
    a = ap.parse_args()

    res = compare(read_vcf(a.calls), read_vcf(a.truth), read_bed(a.bed))
    d = res.as_dict()
    Path(a.json).write_text(json.dumps(d, indent=2))
    with open(a.tsv, "w") as fh:
        fh.write("type\tTP\tFP\tFN\tPPA\tPPV\tGT_mismatch\n")
        for t, m in d["by_type"].items():
            fh.write(f"{t}\t{m['TP']}\t{m['FP']}\t{m['FN']}\t{m['PPA']}\t{m['PPV']}\t{m['GT_mismatch']}\n")
    print(Path(a.tsv).read_text())
    return 0


if __name__ == "__main__":
    sys.exit(main())
