"""Most basic whole-body model (WBM) test, Python side, on Harvey 1.03d.

LP exactly as MATLAB optimizeCbModel builds it for a WBM:
    opt c'v  s.t.  S v (csense) b,  C v (dsense) d,  lb <= v <= ub
Solved with SciPy's HiGHS (dual simplex and interior point, as a cross-check).

Test 1  feasibility: model as shipped (Whole_body_objective_rxn fixed lb=ub=1)
Test 2  max Whole_body_objective_rxn with bounds relaxed to [0, 1000]
Test 3  min Whole_body_objective_rxn with bounds relaxed to [0, 1000]

Usage:  python wbm_basic_python.py path/to/Harvey_1_03d.mat [results.json]
Needs only numpy and scipy (>= 1.9 for HiGHS options used here).
"""
import sys, time, json
import numpy as np
import scipy.io as sio
import scipy.sparse as sp
from scipy.optimize import linprog

path = sys.argv[1]
t0 = time.time()
raw = sio.loadmat(path, squeeze_me=True, struct_as_record=False)
m = raw[[k for k in raw if not k.startswith("__")][0]]
t_load = time.time() - t0

S = sp.csr_matrix(m.S, dtype=float); C = sp.csr_matrix(m.C, dtype=float)
lb = np.asarray(m.lb, float); ub = np.asarray(m.ub, float)
b = np.asarray(m.b, float); d = np.atleast_1d(np.asarray(m.d, float))
def senses(x, n):
    a = np.asarray(x)
    s = "".join(a.ravel().astype(str)) if a.ndim else str(a)
    assert len(s) == n, (len(s), n)
    return np.array(list(s))
csense = senses(m.csense, S.shape[0]); dsense = senses(m.dsense, C.shape[0])
rxns = [str(r) for r in m.rxns]
wbo = rxns.index("Whole_body_objective_rxn")

info = {
    "model": str(getattr(m, "modelID", "Harvey_1_03d")),
    "n_mets": S.shape[0], "n_rxns": S.shape[1], "n_coupling": C.shape[0],
    "nnz_S": int(S.nnz), "nnz_C": int(C.nnz),
    "csense": {k: int(v) for k, v in zip(*np.unique(csense, return_counts=True))},
    "dsense": {k: int(v) for k, v in zip(*np.unique(dsense, return_counts=True))},
    "shipped_wbo_bounds": [lb[wbo], ub[wbo]], "osenseStr": str(m.osenseStr),
    "nonzero_c": int(np.count_nonzero(m.c)), "load_s": round(t_load, 1),
}
print(json.dumps(info, indent=1))

def to_ub_rows(M, rhs, sense):
    """Return (A_eq, b_eq, A_ub, b_ub) for rows with sense E/L/G."""
    E, L, G = sense == "E", sense == "L", sense == "G"
    A_ub = sp.vstack([M[L], -M[G]]).tocsr(); b_ub = np.concatenate([rhs[L], -rhs[G]])
    return M[E], rhs[E], A_ub, b_ub

Se, be, Su, bu = to_ub_rows(S, b, csense)
Ce, de, Cu, du = to_ub_rows(C, d, dsense)
A_eq = sp.vstack([Se, Ce]).tocsr(); b_eq = np.concatenate([be, de])
A_ub = sp.vstack([Su, Cu]).tocsr(); b_ub = np.concatenate([bu, du])

def check(v, lo, hi):
    Sv = S @ v - b
    sviol = np.max(np.abs(Sv[csense == "E"])) if (csense == "E").any() else 0
    Cv = C @ v - d
    cviol = max(0, np.max(np.where(dsense == "L", Cv, np.where(dsense == "G", -Cv, np.abs(Cv)))))
    bviol = max(0, np.max(lo - v), np.max(v - hi))
    return float(sviol), float(cviol), float(bviol)

results = []
def solve(name, cvec, lo, hi, method):
    t0 = time.time()
    r = linprog(cvec, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq,
                bounds=np.column_stack([lo, hi]), method=method)  # default tolerances: 1e-9 made the feasible model "infeasible"
    t = time.time() - t0
    row = {"test": name, "solver": method, "status": r.status, "message": r.message.split(".")[0],
           "time_s": round(t, 1)}
    if r.status == 0:
        row["wbo_flux"] = float(r.x[wbo])
        row["Sv_max_resid"], row["Cv_max_viol"], row["bound_max_viol"] = check(r.x, lo, hi)
        row["n_active_rxns"] = int(np.sum(np.abs(r.x) > 1e-9))
    print(json.dumps(row)); sys.stdout.flush(); results.append(row)
    return r

zero = np.zeros(S.shape[1])
lo2, hi2 = lb.copy(), ub.copy(); lo2[wbo], hi2[wbo] = 0.0, 1000.0
e = np.zeros(S.shape[1]); e[wbo] = 1.0

for method in ("highs-ds", "highs-ipm"):
    solve("1 feasibility (WBO fixed = 1)", zero, lb, ub, method)
    solve("2 max WBO (bounds 0..1000)", -e, lo2, hi2, method)
    solve("3 min WBO (bounds 0..1000)", e, lo2, hi2, method)

out = sys.argv[2] if len(sys.argv) > 2 else "results_wbm_basic_python.json"
with open(out, "w") as fh:
    json.dump({"model": info, "results": results}, fh, indent=1)
print("wrote", out)
