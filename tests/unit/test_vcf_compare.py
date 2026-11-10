"""Unit tests: the evaluation code itself (REQ-011, REQ-012, REQ-009 filter logic).

If the comparison code is wrong, every accuracy claim built on it is wrong,
so the test tool is validated with hand-counted cases.
"""

import gzip

import pytest

from vvdemo.vcf import (Variant, compare, in_regions, merge_intervals,
                        passes_qual_filter, read_bed, read_vcf, variant_type)

HEADER = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"


def vcf(tmp_path, body, name="x.vcf", gz=False):
    p = tmp_path / name
    text = HEADER + body
    if gz:
        with gzip.open(p, "wt") as fh:
            fh.write(text)
    else:
        p.write_text(text)
    return p


@pytest.mark.req("REQ-011")
@pytest.mark.parametrize("ref, alt, expected", [
    ("A", "G", "SNV"), ("A", "AT", "INDEL"), ("AT", "A", "INDEL"), ("AC", "GT", "MNV"),
])
def test_variant_type(ref, alt, expected):
    assert variant_type(ref, alt) == expected


@pytest.mark.req("REQ-011")
def test_read_vcf_splits_multiallelic_and_filters(tmp_path):
    p = vcf(tmp_path,
            "c\t10\t.\tA\tG,T\t50\tPASS\t.\tGT\t1/2\n"
            "c\t20\t.\tC\tT\t5\tLowQual\t.\tGT\t0/1\n"
            "c\t30\t.\tG\t*\t50\tPASS\t.\tGT\t0/1\n", gz=True, name="x.vcf.gz")
    vs = read_vcf(p)
    assert [v.key for v in vs] == [("c", 10, "A", "G"), ("c", 10, "A", "T")]
    assert len(read_vcf(p, pass_only=False)) == 3


@pytest.mark.req("REQ-011")
def test_read_vcf_rejects_malformed_line(tmp_path):
    with pytest.raises(ValueError, match="malformed"):
        read_vcf(vcf(tmp_path, "c\t10\tA\tG\n"))


@pytest.mark.req("REQ-011")
def test_compare_counts_hand_checked_example():
    truth = [Variant("c", 10, "A", "G", "0/1"), Variant("c", 20, "C", "T", "1/1"),
             Variant("c", 30, "G", "GA", "0/1"), Variant("c", 40, "TA", "T", "0/1")]
    calls = [Variant("c", 10, "A", "G", "0/1"),       # TP
             Variant("c", 20, "C", "T", "0/1"),       # TP, wrong genotype
             Variant("c", 30, "G", "GA", "1|0"),      # TP (phased GT normalised)
             Variant("c", 50, "C", "A", "0/1")]       # FP; indel at 40 is FN
    r = compare(calls, truth)
    snv, ind = r.by_type["SNV"], r.by_type["INDEL"]
    assert (snv.tp, snv.fp, snv.fn, snv.gt_mismatch) == (2, 1, 0, 1)
    assert (ind.tp, ind.fp, ind.fn) == (1, 0, 1)
    assert snv.ppa == 1.0 and snv.ppv == pytest.approx(0.6667)
    assert ind.ppa == 0.5 and ind.ppv == 1.0
    assert r.as_dict()["false_negatives"] == [["c", 40, "TA", "T"]]


@pytest.mark.req("REQ-011")
def test_metrics_undefined_when_no_denominator():
    r = compare([], [Variant("c", 1, "A", "G")])
    assert r.by_type["SNV"].ppv is None and r.by_type["SNV"].ppa == 0.0


@pytest.mark.req("REQ-011")
def test_case_insensitive_alleles():
    r = compare([Variant("c", 1, "a", "g")], [Variant("c", 1, "A", "G")])
    assert r.by_type["SNV"].tp == 1


@pytest.mark.req("REQ-012")
@pytest.mark.parametrize("intervals, expected", [
    ([(5, 10), (1, 3), (2, 4)], [(1, 4), (5, 10)]),   # overlap + unsorted
    ([(1, 5), (5, 8)], [(1, 8)]),                    # touching
    ([], []),
])
def test_merge_intervals(intervals, expected):
    assert merge_intervals(intervals) == expected


@pytest.mark.req("REQ-012")
@pytest.mark.parametrize("pos, ref, inside", [
    (101, "A", True),     # first base of interval [100, 200)
    (100, "A", False),    # one base before
    (200, "A", True),     # last base
    (200, "AC", False),   # deletion spans the interval end
    (201, "A", False),
])
def test_in_regions_boundaries(pos, ref, inside):
    assert in_regions(Variant("c", pos, ref, "G"), {"c": [(100, 200)]}) is inside


@pytest.mark.req("REQ-012")
def test_compare_excludes_outside_confident_regions(tmp_path):
    bed = tmp_path / "r.bed"
    bed.write_text("# confident regions\ntrack name=conf\n\nc\t100\t200\n")
    truth = [Variant("c", 150, "A", "G"), Variant("c", 50, "A", "G")]
    calls = [Variant("c", 150, "A", "G"), Variant("c", 300, "C", "T")]
    r = compare(calls, truth, read_bed(bed))
    assert (r.by_type["SNV"].tp, r.by_type["SNV"].fp, r.by_type["SNV"].fn) == (1, 0, 0)


@pytest.mark.req("REQ-009")
@pytest.mark.parametrize("qual, keep", [(19.99, False), (20, True), (20.01, True), (0, False)])
def test_qual_filter_boundaries(qual, keep):
    assert passes_qual_filter(qual, threshold=20) is keep


@pytest.mark.req("REQ-009")
@pytest.mark.parametrize("bad, err", [(None, TypeError), ("30", TypeError), (True, TypeError), (-1, ValueError)])
def test_qual_filter_rejects_invalid(bad, err):
    with pytest.raises(err):
        passes_qual_filter(bad, threshold=20)
