"""Regression tests: accuracy and output must not change unexpectedly.

REQ-004/005 absolute thresholds, baseline drift, REQ-008 critical variants,
REQ-009 output format, REQ-006 repeatability.
"""

import gzip
import hashlib
import json
import subprocess
from pathlib import Path

import pytest
import yaml

from conftest import run_pipeline
from vvdemo.vcf import read_vcf

HERE = Path(__file__).parent
CRITERIA = yaml.safe_load((HERE / "acceptance_criteria.yaml").read_text())
BASELINE = json.loads((HERE / "baseline_metrics.json").read_text())


def norm_stats(stderr: str) -> dict[str, int]:
    """Parse 'Lines total/split/joined/realigned/...: 53/0/0/0/0/0' (field set varies by bcftools version)."""
    for line in stderr.splitlines():
        if line.startswith("Lines") and ":" in line:
            names, values = line.split(None, 1)[1].split(":")
            return dict(zip(names.strip().split("/"), map(int, values.strip().split("/"))))
    raise AssertionError(f"no bcftools norm summary in: {stderr}")


@pytest.fixture(scope="module")
def metrics(baseline_run):
    return json.loads((baseline_run["out"] / "SYNTH1" / "eval" / "metrics.json").read_text())


@pytest.mark.req("REQ-004")
def test_snv_accuracy_meets_threshold(metrics):
    m, t = metrics["by_type"]["SNV"], CRITERIA["thresholds"]["SNV"]
    assert m["PPA"] >= t["min_ppa"], f"SNV PPA {m['PPA']} < {t['min_ppa']}; FN: {metrics['false_negatives']}"
    assert m["PPV"] >= t["min_ppv"], f"SNV PPV {m['PPV']} < {t['min_ppv']}; FP: {metrics['false_positives']}"


@pytest.mark.req("REQ-005")
def test_indel_accuracy_meets_threshold(metrics):
    m, t = metrics["by_type"].get("INDEL"), CRITERIA["thresholds"]["INDEL"]
    assert m is not None, "no indels evaluated at all"
    assert m["PPA"] >= t["min_ppa"], f"INDEL PPA {m['PPA']} < {t['min_ppa']}; FN: {metrics['false_negatives']}"
    assert m["PPV"] is not None and m["PPV"] >= t["min_ppv"], f"INDEL PPV {m['PPV']} < {t['min_ppv']}"


@pytest.mark.req("REQ-004", "REQ-005")
@pytest.mark.parametrize("vtype", ["SNV", "INDEL"])
@pytest.mark.parametrize("metric", ["PPA", "PPV"])
def test_no_drift_from_locked_baseline(metrics, vtype, metric):
    base = BASELINE["by_type"][vtype][metric]
    now = metrics["by_type"].get(vtype, {}).get(metric) or 0.0
    assert base - now <= CRITERIA["max_drop_vs_baseline"], (
        f"{vtype} {metric} dropped from {base} (baseline) to {now}")


@pytest.mark.req("REQ-008")
def test_critical_variants_always_called(baseline_run):
    calls = {v.key: v for v in read_vcf(baseline_run["out"] / "SYNTH1" / "SYNTH1.vcf.gz")}
    missing, wrong_gt = [], []
    for line in (HERE / "critical_variants.tsv").read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        chrom, pos, ref, alt, gt = line.split("\t")[:5]
        key = (chrom, int(pos), ref, alt)
        if key not in calls:
            missing.append(key)
        elif calls[key].gt.replace("|", "/") != gt:
            wrong_gt.append((key, calls[key].gt, gt))
    assert not missing, f"critical variants not called (or not PASS): {missing}"
    assert not wrong_gt, f"critical variants with wrong genotype: {wrong_gt}"


def _records(vcf_gz):
    with gzip.open(vcf_gz, "rt") as fh:
        return [l.rstrip("\n").split("\t") for l in fh if not l.startswith("#")]


@pytest.mark.req("REQ-009")
def test_vcf_sorted_and_biallelic(baseline_run):
    recs = _records(baseline_run["out"] / "SYNTH1" / "SYNTH1.vcf.gz")
    positions = [(r[0], int(r[1])) for r in recs]
    assert positions == sorted(positions), "VCF not sorted"
    assert all("," not in r[4] for r in recs), "multi-allelic record found; expected norm -m-"


@pytest.mark.req("REQ-009")
def test_vcf_already_normalised(baseline_run):
    """Re-normalising must change nothing: proves indels are left-aligned."""
    out, s = baseline_run["out"], "SYNTH1"
    vcf, ref = out / s / f"{s}.vcf.gz", out / "reference" / "ref.fa"
    proc = subprocess.run(["bcftools", "norm", "-f", str(ref), "-c", "e", str(vcf)],
                          capture_output=True, text=True, check=True)
    st = norm_stats(proc.stderr)
    assert st["total"] > 0
    assert (st["split"], st["joined"], st["realigned"]) == (0, 0, 0), proc.stderr


@pytest.mark.req("REQ-009")
def test_low_qual_calls_flagged(baseline_run):
    for r in _records(baseline_run["out"] / "SYNTH1" / "SYNTH1.vcf.gz"):
        qual, filt = float(r[5]), r[6]
        if qual < 20:
            assert filt == "LowQual", r[:7]
        else:
            assert filt == "PASS", r[:7]


@pytest.mark.req("REQ-011")
def test_truth_set_is_itself_normalised(synthetic, tools_available):
    """A truth set in non-canonical representation would make good calls look like FN+FP."""
    proc = subprocess.run(["bcftools", "norm", "-f", str(synthetic["ref"]), "-c", "e", str(synthetic["truth"])],
                          capture_output=True, text=True, check=True)
    st = norm_stats(proc.stderr)
    assert st["total"] == 53, "truth set size changed unexpectedly"
    assert (st["split"], st["joined"], st["realigned"]) == (0, 0, 0), proc.stderr


def _body_md5(vcf_gz):
    with gzip.open(vcf_gz, "rt") as fh:
        body = "".join(l for l in fh if not l.startswith("##"))  # drop headers with dates/commands
    return hashlib.md5(body.encode()).hexdigest()


@pytest.mark.req("REQ-006")
def test_repeatable_across_runs_and_thread_counts(baseline_run, tmp_path):
    ref_md5 = _body_md5(baseline_run["out"] / "SYNTH1" / "SYNTH1.vcf.gz")
    for i, threads in enumerate([1, 4]):
        wd = tmp_path / f"rep{i}"
        proc = run_pipeline(wd, {"threads": threads}, cores=threads)
        assert proc.returncode == 0, proc.stderr[-2000:]
        assert _body_md5(wd / "results" / "SYNTH1" / "SYNTH1.vcf.gz") == ref_md5, (
            f"VCF differs from baseline when run with {threads} thread(s)")
