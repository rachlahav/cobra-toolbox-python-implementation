"""Quick checks for why a whole-body model LP is infeasible.

Usage:  python diagnose_infeasibility.py path/to/Harvey_1_03d.mat

1. Data checks: lb > ub, NaN/Inf, nonzero b/d.
2. Interval-arithmetic check: rows of S v = b and C v (dsense) d that cannot be
   met by ANY v within the bounds, row by row.
3. LP feasibility (HiGHS) with only the 'L' coupling rows, then only the 'G' rows,
   with Whole_body_objective_rxn relaxed to [0, 1000].
"""
import sys, time
import numpy as np
import scipy.io as sio
import scipy.sparse as sp
from scipy.optimize import linprog

raw = sio.loadmat(sys.argv[1], squeeze_me=True, struct_as_record=False)
m = raw[[k for k in raw if not k.startswith("__")][0]]
S = sp.csr_matrix(m.S, dtype=float); C = sp.csr_matrix(m.C, dtype=float)
lb = np.asarray(m.lb, float); ub = np.asarray(m.ub, float)
b = np.asarray(m.b, float); d = np.atleast_1d(np.asarray(m.d, float))
ds = np.array(list("".join(np.asarray(m.dsense).ravel().astype(str))))
rx = [str(r) for r in m.rxns]; mets = [str(x) for x in m.mets]
w = rx.index("Whole_body_objective_rxn")

print("1. data")
print("   lb > ub:", int((lb > ub).sum()))
print("   NaN in lb/ub/d/S/C:", int(np.isnan(lb).sum() + np.isnan(ub).sum() + np.isnan(d).sum()
                                   + np.isnan(S.data).sum() + np.isnan(C.data).sum()))
print("   nonzero b:", int(np.count_nonzero(b)), " nonzero d:", int(np.count_nonzero(d)))
print("   forced flux: lb > 0:", int((lb > 0).sum()), " ub < 0:", int((ub < 0).sum()))
print("   lb range", lb.min(), lb.max(), "| ub range", ub.min(), ub.max())

print("2. rows impossible on their own (interval arithmetic)")
lbf = np.where(np.isfinite(lb), lb, -1e30); ubf = np.where(np.isfinite(ub), ub, 1e30)
def row_range(M):
    pos, neg = M.maximum(0), M.minimum(0)
    return pos @ lbf + neg @ ubf, pos @ ubf + neg @ lbf
lo, hi = row_range(S)
badS = np.flatnonzero((lo > b + 1e-9) | (hi < b - 1e-9))
lo, hi = row_range(C)
badC = np.flatnonzero(((ds == "L") & (lo > d + 1e-9)) | ((ds == "G") & (hi < d - 1e-9)))
print("   S rows:", len(badS), [mets[i] for i in badS[:10]])
print("   C rows:", len(badC), list(badC[:10]))

print("3. LP feasibility with one coupling group at a time (WBO in [0, 1000])")
lb2, ub2 = lb.copy(), ub.copy(); lb2[w], ub2[w] = 0.0, 1000.0
for sense in ("L", "G"):
    A = C[ds == sense] if sense == "L" else -C[ds == sense]
    rhs = d[ds == sense] if sense == "L" else -d[ds == sense]
    t0 = time.time()
    r = linprog(np.zeros(S.shape[1]), A_ub=A, b_ub=rhs, A_eq=S, b_eq=b,
                bounds=np.column_stack([lb2, ub2]), method="highs", options={"time_limit": 300})
    print(f"   only '{sense}' rows ({A.shape[0]}): status {r.status} - {r.message[:60]} ({time.time()-t0:.1f}s)")
