"""Time a Python optimizeCbModel translation on the same cases as
matlab/run_matlab_reference.m, writing the same CSV schema so compare.py can
line the two up.

Usage:  OPTCB_IMPLS=translation python benchmarks/benchmark.py [repeats] [out.csv]
        (names: translation, reference, or any module:function)
"""
from __future__ import annotations

import copy
import csv
import os
import statistics
import sys
import time

import numpy as np
import scipy.io as sio

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "tests"))
from cobra_mat import load_cobra_mat  # noqa: E402
from test_optimize_cb_model import MODELS, load_impl, normalize  # noqa: E402

CASES = [  # (case id, osenseStr, minNorm, weighted g1)
    ("fba_max", "max", 0, False),
    ("l1_one", "max", "one", False),
    ("l1_weighted", "max", "one", True),
]


def main(repeats=10, out="python_results.csv"):
    repeats = int(repeats)
    spec, impl = load_impl()
    rows, fluxes = [], {}
    for name, path in MODELS.items():
        t0 = time.perf_counter()
        base = load_cobra_mat(path)
        load_s = time.perf_counter() - t0
        for case, sense, mn, weighted in CASES:
            m = copy.deepcopy(base)
            if weighted:
                m["g1"] = np.arange(1, m["S"].shape[1] + 1, dtype=float)
            times, res = [], None
            for _ in range(repeats + 1):  # first call is warm-up, not timed
                mm = copy.deepcopy(m)
                t0 = time.perf_counter()
                try:
                    r = impl(mm, sense, mn)
                except NotImplementedError:
                    r = None
                    break
                times.append(time.perf_counter() - t0)
                res = r
            if res is None:
                continue
            times = times[1:]
            s = normalize(res, m["rxns"], m["c"])
            v = np.asarray(s["v"], float)
            fluxes[f"{name}__{case}"] = v
            rows.append({
                "model": name, "case": case, "stat": s["stat"],
                "f": f"{s['f']:.15g}", "sum_abs_v": f"{np.abs(v).sum():.15g}",
                "mean_s": f"{statistics.mean(times):.6f}",
                "sd_s": f"{statistics.stdev(times) if len(times) > 1 else 0:.6f}",
                "median_s": f"{statistics.median(times):.6f}",
                "repeats": len(times), "load_s": f"{load_s:.4f}",
                "impl": spec,
            })
            print(f"{name:10} {case:12} f={s['f']:.10f} sum|v|={np.abs(v).sum():.6f} "
                  f"median {statistics.median(times)*1e3:8.2f} ms")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    sio.savemat(out.replace(".csv", "_fluxes.mat"), fluxes)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main(*sys.argv[1:])
