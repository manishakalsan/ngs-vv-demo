"""Shared fixtures and the requirement-tracing hook.

Every test declares the requirement(s) it verifies with
@pytest.mark.req("REQ-00X"). After the session, results are written to
reports/test_results.json, which scripts/build_traceability.py turns into
the traceability matrix.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import time
from pathlib import Path

import pytest
import yaml

REPO = Path(__file__).resolve().parent.parent
REPORTS = REPO / "reports"
BASE_CONFIG = REPO / "config" / "config.yaml"
TOOLS = ["snakemake", "fastp", "bwa", "samtools", "bcftools"]

_results: list[dict] = []


# ---------- requirement tracing ----------

def pytest_collection_modifyitems(items):
    """Tag each test with its level from the folder it lives in."""
    for item in items:
        level = Path(str(item.fspath)).parent.name
        if level in {"unit", "smoke", "negative", "regression", "performance"}:
            item.add_marker(getattr(pytest.mark, level))


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_makereport(item, call):
    outcome = yield
    rep = outcome.get_result()
    # record the call phase, or a setup failure/skip
    if rep.when == "call" or (rep.when == "setup" and rep.outcome != "passed"):
        reqs = [r for m in item.iter_markers("req") for r in m.args]
        _results.append({
            "test": item.nodeid,
            "level": Path(str(item.fspath)).parent.name,
            "requirements": sorted(set(reqs)),
            "outcome": rep.outcome,
            "duration_s": round(rep.duration, 3),
        })


def pytest_sessionfinish(session, exitstatus):
    if not _results:
        return
    REPORTS.mkdir(exist_ok=True)
    (REPORTS / "test_results.json").write_text(json.dumps({
        "generated": time.strftime("%Y-%m-%d %H:%M:%S %Z"),
        "results": _results,
    }, indent=2))


# ---------- pipeline helpers ----------

def deep_merge(base: dict, over: dict) -> dict:
    out = dict(base)
    for k, v in over.items():
        out[k] = deep_merge(out[k], v) if isinstance(v, dict) and isinstance(out.get(k), dict) else v
    return out


def run_pipeline(workdir: Path, overrides: dict | None = None, cores: int = 2,
                 timeout: int = 600) -> subprocess.CompletedProcess:
    """Run the Snakefile in an isolated work dir with an optional config override."""
    workdir.mkdir(parents=True, exist_ok=True)
    cfg = yaml.safe_load(BASE_CONFIG.read_text())
    # A proposed change (e.g. config/changes/*.yaml) can be layered on top, so the
    # SAME suite verifies the candidate release: VVDEMO_CONFIG_OVERLAY=path pytest
    overlay = os.environ.get("VVDEMO_CONFIG_OVERLAY")
    if overlay:
        cfg = deep_merge(cfg, yaml.safe_load((REPO / overlay).read_text()))
    cfg = deep_merge(cfg, overrides or {})
    cfg["outdir"] = str(workdir / "results")
    cfg_path = workdir / "config.yaml"
    cfg_path.write_text(yaml.safe_dump(cfg))
    return subprocess.run(
        ["snakemake", "-s", str(REPO / "Snakefile"), "--configfile", str(cfg_path),
         "--cores", str(cores), "--directory", str(workdir), "--quiet"],
        capture_output=True, text=True, timeout=timeout,
    )


@pytest.fixture(scope="session")
def tools_available():
    missing = [t for t in TOOLS if shutil.which(t) is None]
    if missing:
        pytest.fail(f"required tools not on PATH: {missing}. Activate the conda env (see README).")


@pytest.fixture(scope="session")
def baseline_run(tmp_path_factory, tools_available):
    """One full pipeline run on the synthetic dataset, shared by smoke and regression tests."""
    wd = tmp_path_factory.mktemp("baseline")
    t0 = time.perf_counter()
    proc = run_pipeline(wd)
    elapsed = time.perf_counter() - t0
    if proc.returncode != 0:
        pytest.fail(f"baseline pipeline run failed:\n{proc.stderr[-3000:]}")
    out = wd / "results"
    REPORTS.mkdir(exist_ok=True)
    shutil.copy(out / "SYNTH1" / "eval" / "metrics.json", REPORTS / "metrics.json")
    shutil.copy(out / "software_versions.txt", REPORTS / "software_versions.txt")
    return {"workdir": wd, "out": out, "sample": "SYNTH1", "elapsed_s": elapsed, "proc": proc}


@pytest.fixture(scope="session")
def synthetic():
    d = REPO / "data" / "synthetic"
    return {"dir": d, "r1": d / "sample_R1.fastq.gz", "r2": d / "sample_R2.fastq.gz",
            "ref": d / "ref.fa", "truth": d / "truth.vcf", "bed": d / "confident.bed"}
