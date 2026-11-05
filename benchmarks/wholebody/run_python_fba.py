"""Basic whole-body FBA parity test, Python side (cobrapy + GLPK).

Builds the same LP as MATLAB's optimizeCbModel:
  max c'v  s.t.  S v (csense) b,  C v (dsense) d,  lb <= v <= ub
then compares the optimum with matlab_result.mat.

Run run_matlab_fba.m first; it writes harvey_lean.mat and matlab_result.mat next
to this script. Usage:  python run_python_fba.py [solver]   (default: glpk)
NOT YET RUN: the MATLAB side has not produced its outputs yet (see
docs/wholebody_basic_test_findings.md).
"""
import sys, time
from pathlib import Path
import numpy as np
import scipy.io as sio
import scipy.sparse as sp
import cobra
from optlang.symbolics import Zero

here = Path(__file__).resolve().parent
m = sio.loadmat(here / "harvey_lean.mat", squeeze_me=True)
ml = sio.loadmat(here / "matlab_result.mat", squeeze_me=True)

S = sp.csc_matrix(m["S"]); C = sp.csr_matrix(m["C"])
lb, ub, c = m["lb"].astype(float), m["ub"].astype(float), m["c"].astype(float)
rxns, mets = [str(r) for r in m["rxns"]], [str(x) for x in m["mets"]]
d, dsense = np.atleast_1d(m["d"]).astype(float), str(m["dsense"])
assert set(str(m["csense"])) == {"E"} and not np.any(m["b"]), "expected S v = 0"

t0 = time.time()
model = cobra.Model("Harvey_1_03d")
model.solver = sys.argv[1] if len(sys.argv) > 1 else "glpk"
model.add_metabolites([cobra.Metabolite(x) for x in mets])
rx = []
for j, rid in enumerate(rxns):
    r = cobra.Reaction(rid, lower_bound=lb[j], upper_bound=ub[j])
    col = S.getcol(j)
    r.add_metabolites({model.metabolites[i]: v for i, v in zip(col.indices, col.data)})
    rx.append(r)
model.add_reactions(rx)
t_build = time.time() - t0

# coupling constraints  C v <= / >= d   (v = forward - reverse)
t0 = time.time()
fwd = [r.forward_variable for r in model.reactions]
rev = [r.reverse_variable for r in model.reactions]
cons = []
for k in range(C.shape[0]):
    lo, hi = (None, d[k]) if dsense[k] == "L" else (d[k], None) if dsense[k] == "G" else (d[k], d[k])
    cons.append(model.problem.Constraint(Zero, lb=lo, ub=hi, name=f"coupling_{k}"))
model.add_cons_vars(cons)
model.solver.update()
for k, con in enumerate(cons):
    row = C.getrow(k)
    coefs = {}
    for j, v in zip(row.indices, row.data):
        coefs[fwd[j]] = v; coefs[rev[j]] = -v
    con.set_linear_coefficients(coefs)
t_coupling = time.time() - t0

model.objective = {model.reactions[int(np.flatnonzero(c)[0])]: 1}
model.objective_direction = "max"
model.solver.configuration.tolerances.feasibility = 1e-9
t0 = time.time(); sol = model.optimize(); t_solve = time.time() - t0

v = sol.fluxes.values
viol_S = np.abs(S @ v).max()
Cv = C @ v
viol_C = np.max(np.where(np.array(list(dsense)) == "L", Cv - d, d - Cv).clip(min=0))
viol_bounds = max((lb - v).clip(min=0).max(), (v - ub).clip(min=0).max())

f_py, f_ml = sol.objective_value, float(ml["f"])
rel = abs(f_py - f_ml) / max(1.0, abs(f_ml))
print(f"model: {len(mets)} mets x {len(rxns)} rxns, {C.shape[0]} coupling rows")
print(f"build {t_build:.0f}s, coupling {t_coupling:.0f}s, solve {t_solve:.1f}s (MATLAB solve {float(ml['time_s']):.1f}s)")
print(f"status  python={sol.status}  matlab stat={int(ml['stat'])}")
print(f"objective  python={f_py:.10g}  matlab={f_ml:.10g}  rel.diff={rel:.2e}")
print(f"python feasibility: |Sv|max={viol_S:.1e}  Cv viol={viol_C:.1e}  bound viol={viol_bounds:.1e}")
vm = np.asarray(ml["v"], float)
print(f"flux vectors: max |v_py - v_ml| = {np.abs(v - vm).max():.3g} (alternate optima allowed; not a pass criterion)")
ok = sol.status == "optimal" and int(ml["stat"]) == 1 and rel < 1e-6
print("RESULT:", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
