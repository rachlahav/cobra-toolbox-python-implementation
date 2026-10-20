"""Tests for a Python translation of COBRA Toolbox optimizeCbModel.

Reference values come from the COBRA Toolbox's own MATLAB test,
test/verifiedTests/analysis/testOptimizeCbModel/testOptimizeCbModel.m,
which runs on iAF1260 with tol = 1e-6.

Every test runs once per implementation in IMPLS:
  translation  pscm.core.optimize_cb_model, via adapters.workbench_translation
  reference    tests/reference/reference_optimize.py (SciPy/HiGHS, no cobrapy):
               a control that shows the test itself is right
Pick a subset with OPTCB_IMPLS=translation (comma-separated names, or any
"module:function"). Implementations are called as impl(model_dict, osenseStr,
minNorm) and may return a dict, an object, or a cobrapy Solution.

Run:  pytest -q tests/                       (or: python tests/test_optimize_cb_model.py)
"""
from __future__ import annotations

import copy
import csv
import importlib
import json
import os
import sys
import time
import traceback

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from cobra_mat import load_cobra_mat  # noqa: E402

TOL = 1e-6  # same tolerance as testOptimizeCbModel.m

# MATLAB reference values (iAF1260, from testOptimizeCbModel.m)
MATLAB_REF = {
    "iAF1260": {
        "fba_f": 0.736700938697735,          # asserted, lines 161/174/189
        "l1_sum_abs": 5.997682160714440e02,  # line 114 (commented out in MATLAB)
        "wl1_sum_abs": 6.003485501480500e02,  # asserted, line 125, g1 = 1:n
    },
}

MODELS = {
    "ecoli_core": os.path.join(HERE, "models", "ecoli_core_model.mat"),
    "iAF1260": os.path.join(HERE, "models", "iAF1260.mat"),
}


# --------------------------------------------------------------------------
# Implementation loading and result normalisation
# --------------------------------------------------------------------------
IMPLS = {
    "translation": "adapters:workbench_translation",
    "reference": "reference.reference_optimize:optimize_cb_model",
}


def selected_impls():
    raw = os.environ.get("OPTCB_IMPLS") or os.environ.get("OPTCB_IMPL") or ",".join(IMPLS)
    return [s.strip() for s in raw.split(",") if s.strip()]


def load_impl(name=None):
    name = name or selected_impls()[0]
    spec = IMPLS.get(name, name)
    mod, fn = spec.split(":")
    return name, getattr(importlib.import_module(mod), fn)


def normalize(res, rxns, c=None):
    """Map dict / attribute object / cobrapy Solution onto MATLAB field names.
    For cobrapy Solutions, f is recomputed as c'v when c is given, because
    e.g. pfba() reports sum|v| as objective_value, not the growth rate."""
    if isinstance(res, dict):
        get = res.get
    else:
        get = lambda k, d=None: getattr(res, k, d)  # noqa: E731
    out = {}
    # cobrapy Solution
    if get("objective_value") is not None and get("f") is None:
        status = get("status")
        out["stat"] = {"optimal": 1, "infeasible": 0, "unbounded": 2}.get(status, -1)
        out["f"] = get("objective_value")
        fl = get("fluxes")
        out["v"] = np.asarray(fl.reindex(rxns).values if hasattr(fl, "reindex") else fl, float)
        sp_, rc = get("shadow_prices"), get("reduced_costs")
        out["y_cobrapy"] = None if sp_ is None else np.asarray(sp_, float)
        out["w_cobrapy"] = None if rc is None else np.asarray(rc, float)
        if c is not None and out["stat"] == 1:
            out["f"] = float(np.asarray(c, float) @ out["v"])
        elif out["stat"] != 1:
            out["f"] = float("nan")
        out["_source"] = "cobrapy.Solution"
        return out
    for k in ("f", "f0", "f1", "f2", "obj", "stat", "origStat", "y", "w", "s", "time"):
        out[k] = get(k)
    v = get("v")
    x = get("x")
    out["v"] = np.asarray(v if v is not None else x, float) if (v is not None or x is not None) else None
    out["x"] = None if x is None else np.asarray(x, float)
    for k in ("y", "w", "s"):
        if out[k] is not None:
            out[k] = np.asarray(out[k], float)
    out["_source"] = "struct-like"
    return out


class Skip(Exception):
    pass


_model_cache: dict = {}


def get_model(name):
    if name not in _model_cache:
        _model_cache[name] = load_cobra_mat(MODELS[name])
    return copy.deepcopy(_model_cache[name])


def run(impl, model, osenseStr="max", minNorm=0):
    try:
        res = impl(model, osenseStr, minNorm)
    except NotImplementedError as e:
        raise Skip(f"not implemented: {e}")
    return normalize(res, model["rxns"], model["c"])


# --------------------------------------------------------------------------
# Tests: each returns a short detail string or raises AssertionError / Skip
# --------------------------------------------------------------------------
def t_fba_status_and_mass_balance(impl, name):
    m = get_model(name)
    s = run(impl, m)
    assert s["stat"] == 1, f"stat={s['stat']}"
    r = np.linalg.norm(m["S"] @ s["v"] - m["b"])
    assert r < TOL, f"||Sv-b||={r:.2e}"
    return f"||Sv-b||={r:.1e}"


def t_fba_objective_matches_matlab(impl, name):
    s = run(impl, get_model(name))
    if name not in MATLAB_REF:
        raise Skip("no MATLAB reference value for this model yet")
    ref = MATLAB_REF[name]["fba_f"]
    d = abs(s["f"] - ref)
    assert d < TOL, f"f={s['f']:.12f} vs MATLAB {ref:.12f} (diff {d:.2e})"
    return f"f={s['f']:.12f}, diff {d:.1e}"


def t_bounds_respected(impl, name):
    m = get_model(name)
    v = run(impl, m)["v"]
    lo = np.max(m["lb"] - v)
    hi = np.max(v - m["ub"])
    assert lo < TOL and hi < TOL, f"lb violation {lo:.2e}, ub violation {hi:.2e}"
    return f"max violation {max(lo, hi, 0):.1e}"


def t_duals_stationarity(impl, name):
    m = get_model(name)
    s = run(impl, m)
    if s.get("y") is None or s.get("w") is None:
        raise Skip("no y/w (COBRA-convention duals) returned")
    osense = -1  # max
    r = np.linalg.norm(osense * m["c"] - m["S"].T @ s["y"] - s["w"], np.inf)
    assert r < 1e-6, f"||osense*c - S'y - w||_inf = {r:.2e} (COBRA sign convention)"
    return f"stationarity residual {r:.1e}"


def t_minimize_sense(impl, name):
    m = get_model(name)
    s = run(impl, m, "min")
    assert s["stat"] == 1, f"stat={s['stat']}"
    # biomass has lb=0 in both models, so min growth is 0
    assert abs(s["f"]) < TOL, f"min growth f={s['f']}"
    return f"f_min={s['f']:.1e}"


def t_l1_preserves_objective(impl, name):
    m = get_model(name)
    f = run(impl, m)["f"]
    s = run(impl, get_model(name), "max", "one")
    assert s["stat"] == 1, f"stat={s['stat']}"
    r = np.linalg.norm(m["S"] @ s["v"] - m["b"])
    assert r < TOL, f"||Sv-b||={r:.2e}"
    d = abs(f - m["c"] @ s["v"])
    assert d < 0.01, f"objective drifted by {d:.2e} (MATLAB allows < 0.01)"
    return f"objective drift {d:.1e}"


def t_l1_total_flux_matches_matlab(impl, name):
    if name not in MATLAB_REF:
        raise Skip("no MATLAB reference value for this model yet")
    s = run(impl, get_model(name), "max", "one")
    tot = np.abs(s["v"]).sum()
    ref = MATLAB_REF[name]["l1_sum_abs"]
    d = abs(tot - ref)
    # MATLAB has this assert commented out; use a looser 1e-5
    assert d < 1e-5, f"sum|v|={tot:.10f} vs MATLAB {ref:.10f} (diff {d:.2e})"
    return f"sum|v|={tot:.8f}, diff {d:.1e}"


def t_weighted_l1_matches_matlab(impl, name):
    if name not in MATLAB_REF:
        raise Skip("no MATLAB reference value for this model yet")
    m = get_model(name)
    m["g1"] = np.arange(1, m["S"].shape[1] + 1, dtype=float)
    s = run(impl, m, "max", "one")
    tot = np.abs(s["v"]).sum()
    ref = MATLAB_REF[name]["wl1_sum_abs"]
    d = abs(tot - ref)
    assert d < TOL, f"sum|v|={tot:.10f} vs MATLAB {ref:.10f} (diff {d:.2e})"
    return f"sum|v|={tot:.8f}, diff {d:.1e}"


def t_f1_is_unweighted_l1(impl, name):
    m = get_model(name)
    m["g1"] = np.arange(1, m["S"].shape[1] + 1, dtype=float)
    s = run(impl, m, "max", "one")
    if s.get("f1") is None:
        raise Skip("no f1 field returned")
    d = abs(s["f1"] - np.abs(s["v"]).sum())
    assert d < TOL, f"f1={s['f1']} but sum|v|={np.abs(s['v']).sum()} (MATLAB f1 is unweighted)"
    return "f1 == sum|v|"


def t_infeasible_returns_stat0_nan(impl, name):
    m = get_model(name)
    j = int(np.flatnonzero(m["c"])[0])
    m["lb"][j] = 1000.0  # force impossible growth
    m["ub"][j] = 1000.0
    s = run(impl, m)
    assert s["stat"] == 0, f"stat={s['stat']} (expected 0 = infeasible)"
    f = s.get("f")
    assert f is None or np.isnan(f), f"f={f} (MATLAB returns NaN)"
    return "stat=0, f=NaN"


def t_input_not_mutated(impl, name):
    m = get_model(name)
    m["g1"] = np.ones(m["S"].shape[1])
    before = copy.deepcopy(m)
    # 'min' differs from the stored osenseStr, so an overwrite is visible
    run(impl, m, "min", "one")
    for k in before:
        a, b = before[k], m.get(k)
        same = (a != b).nnz == 0 if hasattr(a, "nnz") else np.array_equal(np.asarray(a), np.asarray(b))
        assert same, f"model['{k}'] changed (MATLAB passes structs by value)"
    added = set(m) - set(before)
    assert not added, f"new keys added to model: {added}"
    return "model unchanged"


def t_output_fields(impl, name):
    s = run(impl, get_model(name))
    if s["_source"] == "cobrapy.Solution":
        raise Skip("cobrapy Solution: MATLAB field names not expected")
    missing = [k for k in ("f", "v", "x", "y", "w", "s", "stat", "origStat") if s.get(k) is None]
    assert not missing, f"missing fields {missing}"
    assert np.array_equal(s["x"], s["v"]), "x must equal v"
    return "f,v,x,y,w,s,stat,origStat present"


TESTS = [
    t_fba_status_and_mass_balance,
    t_fba_objective_matches_matlab,
    t_bounds_respected,
    t_duals_stationarity,
    t_minimize_sense,
    t_l1_preserves_objective,
    t_l1_total_flux_matches_matlab,
    t_weighted_l1_matches_matlab,
    t_f1_is_unweighted_l1,
    t_infeasible_returns_stat0_nan,
    t_input_not_mutated,
    t_output_fields,
]


# pytest entry points: one test per (implementation, check, model)
def _make_pytest(impl_name, fn, model_name):
    def _t():
        import pytest
        try:
            _, impl = load_impl(impl_name)
        except ImportError as e:
            pytest.skip(f"{impl_name} unavailable: {e}")
        try:
            fn(impl, model_name)
        except Skip as e:
            pytest.skip(str(e))
    _t.__name__ = f"test_{impl_name}__{fn.__name__[2:]}__{model_name}"
    return _t


for _impl in selected_impls():
    _key = _impl.replace(":", "_").replace(".", "_")
    for _fn in TESTS:
        for _name in MODELS:
            globals()[f"test_{_key}__{_fn.__name__[2:]}__{_name}"] = _make_pytest(_impl, _fn, _name)


def main(out_prefix="results"):
    rows = []
    for impl_name in selected_impls():
        _, impl = load_impl(impl_name)
        print(f"\n=== Implementation: {impl_name} ({IMPLS.get(impl_name, impl_name)})\n")
        for name in MODELS:
            for fn in TESTS:
                t0 = time.perf_counter()
                try:
                    detail, status = fn(impl, name), "PASS"
                except Skip as e:
                    detail, status = str(e), "SKIP"
                except AssertionError as e:
                    detail, status = str(e), "FAIL"
                except Exception as e:  # noqa: BLE001
                    detail, status = f"{type(e).__name__}: {e}", "ERROR"
                    traceback.print_exc()
                dt = time.perf_counter() - t0
                rows.append({"impl": impl_name, "model": name, "test": fn.__name__[2:],
                             "status": status, "seconds": round(dt, 4), "detail": detail})
                print(f"{status:5}  {name:10}  {fn.__name__[2:]:32}  {dt:7.3f}s  {detail}")
    counts = {}
    for impl_name in selected_impls():
        counts[impl_name] = {k: sum(r["status"] == k and r["impl"] == impl_name for r in rows)
                             for k in ("PASS", "FAIL", "SKIP", "ERROR")}
        print(f"\n{impl_name:12} " + "  ".join(f"{k}: {v}" for k, v in counts[impl_name].items()))
    with open(f"{out_prefix}.json", "w") as fh:
        json.dump({"counts": counts, "results": rows}, fh, indent=2)
    with open(f"{out_prefix}.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    bad = sum(c["FAIL"] + c["ERROR"] for c in counts.values())
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:]))
