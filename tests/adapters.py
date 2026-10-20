"""Adapters that plug an implementation into the test suite and benchmark.

The suite calls  impl(model_dict, osenseStr, minNorm)  with the output of
cobra_mat.load_cobra_mat (S, b, c, lb, ub, csense, rxns, mets, path and
optionally g1). An adapter turns that into whatever the implementation takes.

    adapters:workbench_translation  pscm.core.optimize_cb_model (cobrapy)
    adapters:cobrapy_builtin        plain cobrapy optimize() / pfba(), as a baseline
    reference.reference_optimize:optimize_cb_model   SciPy/HiGHS port, no cobrapy

Solver: OPTCB_SOLVER=gurobi (etc.) to choose; default GLPK, else SciPy.
"""
from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))  # repo root, for `import pscm`

import copy
from contextlib import contextmanager

_cache: dict = {}
SOLVER = os.environ.get("OPTCB_SOLVER")  # e.g. 'gurobi'; default GLPK, else SciPy


def _maybe_tolerance_shim():
    """OPTCB_TOLERANCE_SHIM=1: let optlang solvers without feasibility/optimality
    tolerances (SciPy, GLPK lacks optimality) accept them as stored, ignored
    values. Sandbox-only workaround so the rest of the translation can be
    exercised; OFF by default so normal runs show the real incompatibility."""
    if os.environ.get("OPTCB_TOLERANCE_SHIM") != "1":
        return
    from optlang.util import SolverTolerances
    if getattr(SolverTolerances, "_pscm_shim", False):
        return
    orig_get, orig_set = SolverTolerances.__getattr__, SolverTolerances.__setattr__

    def _get(self, item):
        try:
            return orig_get(self, item)
        except AttributeError:
            return self.__dict__.setdefault("_shim_values", {}).get(item, 1e-6)

    def _set(self, key, value):
        try:
            orig_set(self, key, value)
        except AttributeError:
            self.__dict__.setdefault("_shim_values", {})[key] = value

    SolverTolerances.__getattr__, SolverTolerances.__setattr__ = _get, _set
    SolverTolerances._pscm_shim = True
    print("NOTE: OPTCB_TOLERANCE_SHIM=1 - missing solver tolerances are faked (sandbox only).")


def _base_cobra_model(path):
    import cobra
    from cobra.io import load_matlab_model

    _maybe_tolerance_shim()
    if SOLVER:
        cobra.Configuration().solver = SOLVER
    elif cobra.Configuration().solver is None or "glpk" not in str(cobra.Configuration().solver):
        try:
            cobra.Configuration().solver = "glpk"
        except Exception:  # GLPK not installed: fall back to SciPy/HiGHS
            cobra.Configuration().solver = "scipy"
    if path not in _cache:
        _cache[path] = load_matlab_model(path)
    return _cache[path]


def _snapshot(cm):
    """Cheap fingerprint of everything optimizeCbModel must leave unchanged."""
    return {
        "bounds": [(r.lower_bound, r.upper_bound) for r in cm.reactions],
        "objective": str(cm.objective.expression),
        "direction": cm.objective.direction,
        "n_variables": len(cm.variables),
        "n_constraints": len(cm.constraints),
        "solver": type(cm.solver).__module__,
        "notes": copy.deepcopy(dict(cm.notes)),
    }


@contextmanager
def cobra_model_for(model_dict):
    """Yield the cached cobrapy Model with this test's bound changes applied.
    Changes are made inside `with cm:` so cobrapy undoes them afterwards; this
    avoids a full model.copy() per call (about a minute for iAF1260 without
    symengine). Bounds are set as (lb, ub) pairs so cobrapy's lb <= ub check
    never trips mid-update. Objective coefficients are not touched (the tests
    don't change c; reaction.objective_coefficient would add to the objective)."""
    cm = _base_cobra_model(model_dict["path"])
    notes_before = copy.deepcopy(dict(cm.notes))
    try:
        with cm:
            for rxn, lb, ub in zip(cm.reactions, model_dict["lb"], model_dict["ub"]):
                if (rxn.lower_bound, rxn.upper_bound) != (lb, ub):
                    rxn.bounds = (float(lb), float(ub))
            yield cm
    finally:
        cm.notes.clear()
        cm.notes.update(notes_before)


def to_cobra_model(model_dict):
    """Independent copy, for callers that keep the model (slow for big models)."""
    with cobra_model_for(model_dict) as cm:
        return cm.copy()


def workbench_translation(model_dict, osenseStr="max", minNorm=0):
    """pscm.core.optimize_cb_model on a cobrapy model. Raises AssertionError if
    the call leaves the model changed (MATLAB passes structs by value)."""
    from pscm.core.optimize_cb_model import MATLAB_FIELDS_KEY, optimize_cb_model

    with cobra_model_for(model_dict) as cm:
        if "g1" in model_dict:  # the translation reads MATLAB-only fields from model.notes
            cm.notes.setdefault(MATLAB_FIELDS_KEY, {})["g1"] = np.asarray(model_dict["g1"], float)
        before = _snapshot(cm)
        sol = optimize_cb_model(cm, osenseStr, minNorm)
        after = _snapshot(cm)
        changed = [k for k in before if repr(before[k]) != repr(after[k])]
        assert not changed, f"optimize_cb_model changed the cobrapy model: {changed}"
        return sol


def cobrapy_builtin(model_dict, osenseStr="max", minNorm=0):
    """Plain cobrapy as a baseline: model.optimize() for FBA, pfba for 'one'.
    pfba minimises unweighted sum|v| only, so the weighted test is expected to fail."""
    from cobra.flux_analysis import pfba

    with cobra_model_for(model_dict) as cm:
        cm.objective_direction = osenseStr
        return pfba(cm) if minNorm == "one" else cm.optimize()
