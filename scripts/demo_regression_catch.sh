#!/usr/bin/env bash
# Demonstrates the regression suite catching a harmful change (CR-001).
#
#   1. run regression tests on the released config  -> must PASS
#   2. run the SAME tests with the change applied   -> must FAIL
#   3. write a before/after table to reports/CR-001_impact.md
#
# Exit code 0 means "the suite caught the change", which is the desired outcome.
set -uo pipefail
cd "$(dirname "$0")/.."
CHANGE=config/changes/v1.1_skip_indels.yaml
mkdir -p reports

echo "== Baseline (v1.0.0) =="
pytest tests/smoke tests/regression -q -p no:cacheprovider
base_rc=$?
cp reports/metrics.json reports/metrics_v1.0.json

echo
echo "== Candidate with $CHANGE =="
VVDEMO_CONFIG_OVERLAY=$CHANGE pytest tests/smoke tests/regression -q -p no:cacheprovider \
  | tee reports/CR-001_pytest.txt
cand_rc=${PIPESTATUS[0]}
cp reports/metrics.json reports/metrics_v1.1-rc1.json
# restore baseline test results for the traceability matrix
pytest tests -q -p no:cacheprovider >/dev/null

python - <<'PY'
import json
a = json.load(open("reports/metrics_v1.0.json"))["by_type"]
b = json.load(open("reports/metrics_v1.1-rc1.json"))["by_type"]
rows = ["| Class | Metric | v1.0.0 (released) | v1.1.0-rc1 (CR-001) | Change |",
        "| --- | --- | --: | --: | --: |"]
for t in ("SNV", "INDEL"):
    for m in ("TP", "FN", "FP", "PPA", "PPV"):
        x, y = a.get(t, {}).get(m), b.get(t, {}).get(m)
        d = "" if x is None or y is None else f"{y - x:+.4g}"
        rows.append(f"| {t} | {m} | {x} | {y} | {d} |")
failed = [l.split(" - ")[0].replace("FAILED ", "").strip() for l in open("reports/CR-001_pytest.txt") if l.startswith("FAILED")]
out = ["# CR-001 impact assessment: skip indel calling", "",
       "Change: `caller.call_indels: false` (proposed to speed up `bcftools mpileup`).", "",
       *rows, "", f"Regression tests failed: {len(failed)}", "",
       *[f"- `{f}`" for f in failed], "",
       "**Decision: REJECT.** REQ-005 (indel accuracy) and REQ-008 (critical variants) are violated."]
open("reports/CR-001_impact.md", "w").write("\n".join(out) + "\n")
print("\n".join(out))
PY

if [ "$base_rc" -eq 0 ] && [ "$cand_rc" -ne 0 ]; then
  echo; echo "RESULT: regression suite CAUGHT the harmful change (expected)."; exit 0
fi
echo; echo "RESULT: unexpected (baseline rc=$base_rc, candidate rc=$cand_rc)"; exit 1
