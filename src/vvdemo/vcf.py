"""Minimal VCF reading, filtering and truth comparison (REQ-004, REQ-005).

Comparison is by normalised representation (CHROM, POS, REF, ALT) after the
pipeline runs `bcftools norm -m- -f ref`, restricted to confident regions.

Limitation (documented in docs/soup_list.md and the verification report):
this is not haplotype-aware. Complex variants written differently can be
miscounted. For real GIAB benchmarking use hap.py / rtg vcfeval; see
scripts/fetch_giab.sh and docs/giab_mode.md.
"""

from __future__ import annotations

import gzip
from dataclasses import dataclass, field
from pathlib import Path


@dataclass(frozen=True)
class Variant:
    chrom: str
    pos: int  # 1-based
    ref: str
    alt: str
    gt: str = "./."

    @property
    def key(self) -> tuple[str, int, str, str]:
        return (self.chrom, self.pos, self.ref.upper(), self.alt.upper())

    @property
    def vtype(self) -> str:
        return variant_type(self.ref, self.alt)


def variant_type(ref: str, alt: str) -> str:
    if len(ref) == 1 and len(alt) == 1:
        return "SNV"
    if len(ref) != len(alt):
        return "INDEL"
    return "MNV"


def _open(path):
    path = str(path)
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def read_vcf(path: str | Path, pass_only: bool = True) -> list[Variant]:
    """Read a VCF into Variant records, splitting multi-allelic ALTs."""
    out: list[Variant] = []
    with _open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                raise ValueError(f"{path}: malformed VCF line with {len(f)} columns")
            chrom, pos, _id, ref, alts, _qual, filt = f[:7]
            if pass_only and filt not in ("PASS", "."):
                continue
            gt = "./."
            if len(f) >= 10:
                fmt = f[8].split(":")
                vals = f[9].split(":")
                if "GT" in fmt:
                    gt = vals[fmt.index("GT")]
            for alt in alts.split(","):
                if alt in (".", "*"):
                    continue
                out.append(Variant(chrom, int(pos), ref, alt, gt))
    return out


def read_bed(path: str | Path) -> dict[str, list[tuple[int, int]]]:
    """BED intervals (0-based half-open), merged and sorted per contig."""
    raw: dict[str, list[tuple[int, int]]] = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            c, s, e = line.split("\t")[:3]
            raw.setdefault(c, []).append((int(s), int(e)))
    return {c: merge_intervals(iv) for c, iv in raw.items()}


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Merge overlapping or touching half-open intervals."""
    if not intervals:
        return []
    ivs = sorted(intervals)
    merged = [ivs[0]]
    for s, e in ivs[1:]:
        ps, pe = merged[-1]
        if s <= pe:
            merged[-1] = (ps, max(pe, e))
        else:
            merged.append((s, e))
    return merged


def in_regions(v: Variant, regions: dict[str, list[tuple[int, int]]]) -> bool:
    """True if the variant's REF span lies fully inside one confident interval."""
    start0 = v.pos - 1
    end0 = start0 + len(v.ref)
    for s, e in regions.get(v.chrom, []):
        if s <= start0 and end0 <= e:
            return True
        if s > start0:
            break
    return False


def passes_qual_filter(qual: float | None, threshold: float) -> bool:
    """QUAL filter used by the pipeline: keep calls with QUAL >= threshold.

    Raises ValueError for missing/negative QUAL so a broken caller output
    can't be silently filtered to nothing.
    """
    if qual is None or isinstance(qual, bool) or not isinstance(qual, (int, float)):
        raise TypeError(f"QUAL must be a number, got {qual!r}")
    if qual < 0 or threshold < 0:
        raise ValueError("QUAL and threshold must be >= 0")
    return qual >= threshold


@dataclass
class ClassMetrics:
    tp: int = 0
    fp: int = 0
    fn: int = 0
    gt_mismatch: int = 0

    @property
    def ppa(self) -> float | None:  # sensitivity / recall
        d = self.tp + self.fn
        return round(self.tp / d, 4) if d else None

    @property
    def ppv(self) -> float | None:  # precision
        d = self.tp + self.fp
        return round(self.tp / d, 4) if d else None

    def as_dict(self) -> dict:
        return {"TP": self.tp, "FP": self.fp, "FN": self.fn,
                "GT_mismatch": self.gt_mismatch, "PPA": self.ppa, "PPV": self.ppv}


@dataclass
class Comparison:
    by_type: dict[str, ClassMetrics] = field(default_factory=dict)
    fp_calls: list[Variant] = field(default_factory=list)
    fn_calls: list[Variant] = field(default_factory=list)

    def as_dict(self) -> dict:
        return {
            "by_type": {k: v.as_dict() for k, v in sorted(self.by_type.items())},
            "false_positives": [list(v.key) for v in self.fp_calls],
            "false_negatives": [list(v.key) for v in self.fn_calls],
        }


def _norm_gt(gt: str) -> str:
    alleles = gt.replace("|", "/").split("/")
    return "/".join(sorted(alleles))


def compare(calls: list[Variant], truth: list[Variant],
            regions: dict[str, list[tuple[int, int]]] | None = None) -> Comparison:
    """Count TP/FP/FN per variant class inside confident regions."""
    if regions is not None:
        calls = [v for v in calls if in_regions(v, regions)]
        truth = [v for v in truth if in_regions(v, regions)]
    call_map = {v.key: v for v in calls}
    truth_map = {v.key: v for v in truth}
    res = Comparison()
    for key, t in truth_map.items():
        m = res.by_type.setdefault(t.vtype, ClassMetrics())
        if key in call_map:
            m.tp += 1
            if _norm_gt(call_map[key].gt) != _norm_gt(t.gt):
                m.gt_mismatch += 1
        else:
            m.fn += 1
            res.fn_calls.append(t)
    for key, c in call_map.items():
        if key not in truth_map:
            res.by_type.setdefault(c.vtype, ClassMetrics()).fp += 1
            res.fp_calls.append(c)
    return res
