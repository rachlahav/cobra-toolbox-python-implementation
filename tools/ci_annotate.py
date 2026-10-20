"""Turn a results JSON from tests/test_optimize_cb_model.py into GitHub
Actions annotations, so results can be read from the API (check-run
annotations) as well as the web UI.

Usage: python tools/ci_annotate.py results.json "label"

GitHub keeps at most 10 annotations of each level per step, so failures are
grouped by message: one ::error per distinct problem, listing the tests.
"""
import json
import sys
from collections import defaultdict


def esc(s):
    return s.replace("%", "%25").replace("\r", "%0D").replace("\n", "%0A")


def main(path, label):
    data = json.load(open(path))
    for impl, c in data["counts"].items():
        print(f"::notice title={label} {impl}::" + esc(
            f"PASS {c['PASS']}  FAIL {c['FAIL']}  SKIP {c['SKIP']}  ERROR {c['ERROR']}"))
    groups = defaultdict(list)
    for r in data["results"]:
        if r["status"] in ("FAIL", "ERROR"):
            groups[(r["impl"], r["status"], r["detail"][:300])].append(f"{r['model']}/{r['test']}")
    for i, ((impl, status, detail), tests) in enumerate(sorted(groups.items(), key=lambda kv: -len(kv[1]))):
        if i == 9:
            print(f"::error title={label}::" + esc(f"{len(groups) - 9} more distinct failures not shown"))
            break
        print(f"::error title={label} {impl} {status} x{len(tests)}::" + esc(f"{detail}\nTests: {', '.join(tests)}"))
    passed = [f"{r['model']}/{r['test']}" for r in data["results"] if r["status"] == "PASS"]
    if passed:
        print(f"::notice title={label} passed::" + esc(", ".join(passed)[:3000]))


if __name__ == "__main__":
    main(*sys.argv[1:3])
