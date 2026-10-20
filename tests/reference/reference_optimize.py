"""Reference Python port of COBRA Toolbox optimizeCbModel (LP modes only).

Solver: HiGHS through scipy.optimize.linprog (no cobrapy/optlang needed).
Covers minNorm in {0/None, 'one'} with optional model['g1'] weights, which is
what testOptimizeCbModel.m checks with fixed reference numbers. QP and
cardinality modes ('zero', 'optimizeCardinality', numeric minNorm) are not
implemented here and raise NotImplementedError.

Return value mirrors the MATLAB solution struct field names:
    f, f0, f1, f2, obj, v, x, y, w, s, stat, origStat, time, solver
Stat codes follow COBRA: 1 optimal, 0 infeasible, 2 unbounded, -1 other.
"""
from __future__ import annotations

import time

import numpy as np
import scipy.sparse as sp
from scipy.optimize import linprog

_HIGHS_STAT = {0: 1, 2: 0, 3: 2}  # linprog status -> COBRA stat (else -1)


def _split_constraints(A, b, csense):
    """Turn rows with csense E/L/G into linprog's A_eq / A_ub form."""
    csense = np.array(list(csense))
    eq, le, ge = csense == "E", csense == "L", csense == "G"
    A = sp.csr_matrix(A)
    A_eq, b_eq = A[eq], b[eq]
    A_ub = sp.vstack([A[le], -A[ge]]).tocsr()
    b_ub = np.concatenate([b[le], -b[ge]])
    return A_eq, b_eq, A_ub, b_ub, eq, le, ge


def _solve(c, A, b, csense, lb, ub, osense):
    """Solve min/max c'x s.t. A x (E/L/G) b, lb<=x<=ub. osense: -1 max, 1 min."""
    A_eq, b_eq, A_ub, b_ub, eq, le, ge = _split_constraints(A, b, csense)
    res = linprog(
        osense * c,
        A_ub=A_ub if A_ub.shape[0] else None,
        b_ub=b_ub if A_ub.shape[0] else None,
        A_eq=A_eq if A_eq.shape[0] else None,
        b_eq=b_eq if A_eq.shape[0] else None,
        bounds=np.column_stack([lb, ub]),
        method="highs",
        options={"presolve": True, "primal_feasibility_tolerance": 1e-9,
                 "dual_feasibility_tolerance": 1e-9},
    )
    stat = _HIGHS_STAT.get(res.status, -1)
    out = {"stat": stat, "origStat": res.status, "message": res.message}
    if stat != 1:
        return out
    x = res.x
    # Duals in COBRA convention (solveCobraLP.m): osense*c - A'*y - w = 0.
    # HiGHS solves min(osense*c), so its marginals already satisfy this.
    m = A.shape[0]
    y = np.zeros(m)
    if eq.any():
        y[eq] = res.eqlin.marginals
    if A_ub.shape[0]:
        mu = res.ineqlin.marginals
        nle = int(le.sum())
        y[le] = mu[:nle]
        y[ge] = -mu[nle:]
    w = res.lower.marginals + res.upper.marginals
    out.update(x=x, obj=float(c @ x), y=y, w=w,
               slack=b - A @ x)
    return out


def optimize_cb_model(model: dict, osenseStr: str = "max", minNorm=0,
                      allowLoops: bool = True, param: dict | None = None):
    if not allowLoops:
        raise NotImplementedError("allowLoops=False (loopless MILP) not ported")
    if minNorm in (0, None, False):
        minNorm = None
    if minNorm is not None and minNorm != "one":
        raise NotImplementedError(f"minNorm={minNorm!r} not ported")

    t0 = time.perf_counter()
    S, b, c = model["S"], model["b"], model["c"]
    lb, ub, csense = model["lb"], model["ub"], model["csense"]
    nMets, nRxns = S.shape
    osense = -1 if (osenseStr or model.get("osenseStr", "max")) == "max" else 1

    sol = _solve(c, S, b, csense, lb, ub, osense)
    objectiveLP = sol.get("obj")

    if sol["stat"] == 1 and minNorm == "one":
        # Mirrors optimizeCbModel.m lines 643-674 exactly.
        I = sp.identity(nRxns, format="csr")
        Z = sp.csr_matrix((nMets, 2 * nRxns))
        A2 = sp.vstack([
            sp.hstack([S, Z]),
            sp.hstack([I, I, sp.csr_matrix((nRxns, nRxns))]),   # v + vf >= 0
            sp.hstack([-I, sp.csr_matrix((nRxns, nRxns)), I]),  # -v + vr >= 0
            sp.hstack([sp.csr_matrix(c), sp.csr_matrix((1, 2 * nRxns))]),
        ]).tocsr()
        g1 = np.asarray(model["g1"], float).ravel() if "g1" in model else np.ones(nRxns)
        c2 = np.concatenate([np.zeros(nRxns), g1, g1])
        lb2 = np.concatenate([lb, np.zeros(2 * nRxns)])
        ub2 = np.concatenate([ub, np.full(2 * nRxns, np.inf)])
        b2 = np.concatenate([b, np.zeros(2 * nRxns), [objectiveLP]])
        cs2 = csense + "G" * (2 * nRxns) + ("G" if osense == -1 else "L")
        sol2 = _solve(c2, A2, b2, cs2, lb2, ub2, 1)
        sol2["x"] = sol2.get("x")
        sol = sol2

    elapsed = time.perf_counter() - t0
    res = {"stat": sol["stat"], "origStat": sol["origStat"], "time": elapsed,
           "solver": "highs (scipy)", "f0": np.nan, "f1": np.nan, "f2": np.nan}
    if sol["stat"] in (1, 3):
        x = sol["x"]
        v = x[:nRxns]
        res["f"] = float(c @ v)
        res["f1"] = float(np.abs(v).sum()) if minNorm == "one" else float(c @ v)
        res["obj"] = np.nansum([res["f0"], res["f1"], res["f2"]])
        res["v"] = res["x"] = v
        res["y"] = sol["y"][:nMets]
        res["w"] = sol["w"][:nRxns]
        res["s"] = sol["slack"][:nMets]
    else:
        res.update(f=np.nan, obj=np.nan, v=np.full(nRxns, np.nan))
        res["x"] = res["v"]
    return res
