"""Performance test (REQ-007). Records runtime to reports/performance.json.

The threshold is generous because CI runners vary; the recorded number is
what matters for trend tracking across releases.
"""

import json

import pytest

from conftest import REPORTS

MAX_WALL_S = 120


@pytest.mark.req("REQ-007")
def test_runtime_within_limit(baseline_run):
    elapsed = baseline_run["elapsed_s"]
    REPORTS.mkdir(exist_ok=True)
    (REPORTS / "performance.json").write_text(json.dumps(
        {"dataset": "synthetic 30 kb, 6,000 read pairs", "cores": 2,
         "wall_time_s": round(elapsed, 2), "limit_s": MAX_WALL_S}, indent=2))
    assert elapsed <= MAX_WALL_S, f"pipeline took {elapsed:.1f}s > {MAX_WALL_S}s"
