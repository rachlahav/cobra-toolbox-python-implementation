"""Compare MATLAB and Python runs produced by run_matlab_reference.m and
benchmark.py.

Usage: python compare.py matlab/matlab_results.csv python_results.csv
"""
from __future__ import annotations

import csv
import os
import sys

import numpy as np
import scipy.io as sio

TOL = 1e-6


def read(path):
    with open(path) as fh:
        return {(r["model"], r["case"]): r for r in csv.DictReader(fh)}


def fluxes(csv_path):
    p = csv_path.replace(".csv", "_fluxes.mat")
    if not os.path.exists(p):
        return {}
    d = sio.loadmat(p, squeeze_me=True)
    return {k: np.asarray(v, float).ravel() for k, v in d.items() if not k.startswith("__")}


def main(mat_csv, py_csv):
    M, P = read(mat_csv), read(py_csv)
    Mf, Pf = fluxes(mat_csv), fluxes(py_csv)
    print(f"{'model':10} {'case':12} {'|df|':>9} {'|d sum|v||':>11} {'max|dv|':>9} "
          f"{'MATLAB ms':>10} {'Python ms':>10} {'speedup':>8}  verdict")
    bad = 0
    for key in sorted(set(M) | set(P)):
        if key not in M or key not in P:
            print(f"{key[0]:10} {key[1]:12} missing on {'MATLAB' if key not in M else 'Python'} side")
            continue
        m, p = M[key], P[key]
        df = abs(float(m["f"]) - float(p["f"]))
        ds = abs(float(m["sum_abs_v"]) - float(p["sum_abs_v"]))
        fk = f"{key[0]}__{key[1]}"
        dv = np.max(np.abs(Mf[fk] - Pf[fk])) if fk in Mf and fk in Pf else float("nan")
        tm, tp = float(m["median_s"]) * 1e3, float(p["median_s"]) * 1e3
        ok = m["stat"] == p["stat"] and df < TOL and (key[1] == "fba_max" or ds < TOL)
        # Flux vectors may legitimately differ for plain FBA (alternate optima),
        # so max|dv| is informational there.
        verdict = "OK" if ok else "MISMATCH"
        bad += not ok
        print(f"{key[0]:10} {key[1]:12} {df:9.1e} {ds:11.1e} {dv:9.1e} "
              f"{tm:10.2f} {tp:10.2f} {tm / tp:7.2f}x  {verdict}")
    print("\nspeedup > 1 means Python was faster (median of timed repeats).")
    print("Plain FBA flux vectors can differ at equal objective (alternate optima);"
          " compare f there, and sum|v| for the L1 cases.")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:3]))
