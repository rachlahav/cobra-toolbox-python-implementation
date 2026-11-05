"""Most basic whole-body model (WBM) test: FBA on Harvey 1.03d.

The LP is the one MATLAB optimizeCbModel builds for a WBM:
    opt c'v  s.t.  S v = b,  C v (dsense) d,  lb <= v <= ub
solved with SciPy's HiGHS (no cobrapy needed).

Harvey is not stored in this repo; it ships with the COBRA Toolbox at
test/models/mat/Harvey_1_03d.mat. Point HARVEY_MAT at it, or keep the toolbox at
~/cobratoolbox. Without the file every test here is skipped (e.g. on CI).

    HARVEY_MAT=/path/to/Harvey_1_03d.mat pytest -q tests/test_wholebody_basic.py

Findings (2026-11-05, docs/wholebody_basic_test_findings.md): the model as
shipped IS feasible with HiGHS's default settings. With Whole_body_objective_rxn
relaxed to [0, 1000], HiGHS presolve wrongly reports "infeasible" in < 1 s, even
though the as-shipped solution also satisfies the relaxed problem (checked below
without a solver). Those solves are marked xfail (known solver issue); if one
starts passing, pytest reports XPASS: update the findings doc and the marker.
"""
from __future__ import annotations

import os

import numpy as np
import pytest
import scipy.sparse as sp

sio = pytest.importorskip("scipy.io")
linprog = pytest.importorskip("scipy.optimize").linprog

HARVEY = os.environ.get(
    "HARVEY_MAT",
    os.path.expanduser("~/cobratoolbox/test/models/mat/Harvey_1_03d.mat"),
)
WBO = "Whole_body_objective_rxn"
FEAS_TOL = 1e-6

PRESOLVE_XFAIL = pytest.mark.xfail(
    reason="HiGHS presolve falsely reports the relaxed Harvey LP infeasible "
           "(see docs/wholebody_basic_test_findings.md)",
    strict=False,
)

pytestmark = pytest.mark.skipif(
    not os.path.isfile(HARVEY),
    reason=f"Harvey model not found at {HARVEY} (set HARVEY_MAT)",
)


def _senses(x, n):
    a = np.asarray(x)
    s = "".join(a.ravel().astype(str)) if a.ndim else str(a)
    assert len(s) == n, f"sense vector has {len(s)} entries, expected {n}"
    return np.array(list(s))


@pytest.fixture(scope="module")
def harvey():
    raw = sio.loadmat(HARVEY, squeeze_me=True, struct_as_record=False)
    m = raw[[k for k in raw if not k.startswith("__")][0]]
    S = sp.csr_matrix(m.S, dtype=float)
    C = sp.csr_matrix(m.C, dtype=float)
    rxns = [str(r) for r in m.rxns]
    return {
        "S": S, "C": C,
        "lb": np.asarray(m.lb, float), "ub": np.asarray(m.ub, float),
        "b": np.asarray(m.b, float), "d": np.atleast_1d(np.asarray(m.d, float)),
        "c": np.asarray(m.c, float),
        "csense": _senses(m.csense, S.shape[0]),
        "dsense": _senses(m.dsense, C.shape[0]),
        "rxns": rxns, "wbo": rxns.index(WBO),
    }


def _solve(h, cvec, lb, ub):
    S, C, d, ds = h["S"], h["C"], h["d"], h["dsense"]
    assert set(h["csense"]) == {"E"}
    A_ub = sp.vstack([C[ds == "L"], -C[ds == "G"]]).tocsr()
    b_ub = np.concatenate([d[ds == "L"], -d[ds == "G"]])
    A_eq = C[ds == "E"] if (ds == "E").any() else None
    A_eq = S if A_eq is None else sp.vstack([S, A_eq]).tocsr()
    b_eq = np.concatenate([h["b"], d[ds == "E"]])
    return linprog(cvec, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq,
                   bounds=np.column_stack([lb, ub]), method="highs")


def _assert_feasible(h, v, lb, ub):
    S, C, d, ds = h["S"], h["C"], h["d"], h["dsense"]
    assert np.abs(S @ v - h["b"]).max() < FEAS_TOL, "mass balance violated"
    Cv = C @ v - d
    viol = np.where(ds == "L", Cv, np.where(ds == "G", -Cv, np.abs(Cv)))
    assert viol.max() < FEAS_TOL, "coupling constraint violated"
    assert (lb - v).max() < FEAS_TOL and (v - ub).max() < FEAS_TOL, "bounds violated"


def _relaxed(h):
    lb, ub = h["lb"].copy(), h["ub"].copy()
    lb[h["wbo"]], ub[h["wbo"]] = 0.0, 1000.0
    return lb, ub


# --- model loads as MATLAB stores it -----------------------------------------

def test_dimensions(harvey):
    assert harvey["S"].shape == (56452, 81094)
    assert harvey["C"].shape == (103600, 81094)
    assert harvey["S"].nnz == 283529
    assert harvey["C"].nnz == 209602


def test_constraint_senses(harvey):
    assert set(harvey["csense"]) == {"E"}
    assert int((harvey["dsense"] == "L").sum()) == 70009
    assert int((harvey["dsense"] == "G").sum()) == 33591
    assert not np.any(harvey["b"]) and not np.any(harvey["d"])


def test_whole_body_objective_as_shipped(harvey):
    w = harvey["wbo"]
    assert harvey["lb"][w] == 1 and harvey["ub"][w] == 1
    assert not np.any(harvey["c"]), "shipped model has no objective set"


def test_data_consistent(harvey):
    lb, ub = harvey["lb"], harvey["ub"]
    assert not np.any(lb > ub)
    for k in ("lb", "ub", "b", "d"):
        assert not np.isnan(harvey[k]).any(), f"NaN in {k}"
    assert not np.isnan(harvey["S"].data).any() and not np.isnan(harvey["C"].data).any()


# --- FBA ---------------------------------------------------------------------

@pytest.fixture(scope="module")
def shipped_solution(harvey):
    """Feasibility LP on the model as shipped (WBO fixed at 1). ~3 min."""
    return _solve(harvey, np.zeros(harvey["S"].shape[1]), harvey["lb"], harvey["ub"])


def test_feasible_as_shipped(harvey, shipped_solution):
    """Whole_body_objective_rxn fixed at 1: the model admits a flux."""
    r = shipped_solution
    assert r.status == 0, r.message
    _assert_feasible(harvey, r.x, harvey["lb"], harvey["ub"])
    assert abs(r.x[harvey["wbo"]] - 1) < FEAS_TOL


def test_relaxed_problem_is_feasible(harvey, shipped_solution):
    """No solver: the as-shipped flux also satisfies WBO in [0, 1000], so the
    relaxed problem is feasible and any 'infeasible' answer for it is wrong."""
    assert shipped_solution.status == 0, shipped_solution.message
    lb, ub = _relaxed(harvey)
    _assert_feasible(harvey, shipped_solution.x, lb, ub)


@PRESOLVE_XFAIL
def test_max_whole_body_objective(harvey):
    """Max Whole_body_objective_rxn with bounds [0, 1000]: optimum >= 1."""
    lb, ub = _relaxed(harvey)
    e = np.zeros(harvey["S"].shape[1]); e[harvey["wbo"]] = 1.0
    r = _solve(harvey, -e, lb, ub)
    assert r.status == 0, r.message
    _assert_feasible(harvey, r.x, lb, ub)
    assert r.x[harvey["wbo"]] >= 1 - FEAS_TOL


@PRESOLVE_XFAIL
def test_min_whole_body_objective(harvey):
    """Min Whole_body_objective_rxn with bounds [0, 1000]: optimum <= max."""
    lb, ub = _relaxed(harvey)
    e = np.zeros(harvey["S"].shape[1]); e[harvey["wbo"]] = 1.0
    r = _solve(harvey, e, lb, ub)
    assert r.status == 0, r.message
    _assert_feasible(harvey, r.x, lb, ub)
