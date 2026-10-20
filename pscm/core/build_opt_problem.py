# Exported from the PSCM workbench (record 'buildOptProblemFromModel', status: Drafting, updated 2026-10-05T17:45:03.393Z).
# Edit in the workbench, then re-run tools/export_workbench.py.
"""Build the MATLAB-style LP from a cobrapy whole-body model.

MATLAB source: src/base/solvers/buildOptProblemFromModel.m (cobratoolbox),
commit: TODO (record the hash the port is checked against).
The .m file is not in the synced project files; this port follows the LP
contract documented in src/analysis/FBA/optimizeCbModel.m:

    max/min  c' v
    s.t.     S v  (csense)  b      one row per metabolite
             C v  (dsense)  d      one row per coupling constraint
             lb <= v <= ub

cobrapy keeps the LP inside the solver

Status: Draft. 
"""
from dataclasses import dataclass, field
from typing import Dict, List, Sequence, Tuple

import numpy as np
import scipy.sparse as sp

# Name of the model.notes key listing coupling constraint ids (MATLAB model.ctrs),
# in MATLAB row order. Set by pscm.core.coupling.add_coupling_constraints.
COUPLING_KEY = "coupling_constraint_ids"


@dataclass
class LPProblem:
    """Mirror of the MATLAB LPproblem struct (rows = mets then ctrs, cols = rxns)."""
    A: sp.csr_matrix
    b: np.ndarray
    c: np.ndarray
    lb: np.ndarray
    ub: np.ndarray
    csense: np.ndarray            # array of 'E' / 'L' / 'G'
    osense: int                   # -1 maximise, +1 minimise (MATLAB convention)
    rxn_ids: List[str] = field(default_factory=list)
    met_ids: List[str] = field(default_factory=list)
    ctr_ids: List[str] = field(default_factory=list)

    @property
    def n_mets(self) -> int:
        return len(self.met_ids)

    @property
    def C(self) -> sp.csr_matrix:
        """The coupling block, equal to MATLAB model.C."""
        return self.A[self.n_mets:, :]


def _row_sense(lb, ub, name: str) -> Tuple[str, float]:
    """Translate an optlang constraint's (lb, ub) into MATLAB (sense, rhs)."""
    lb = None if lb is None or np.isneginf(lb) else float(lb)
    ub = None if ub is None or np.isposinf(ub) else float(ub)
    if lb is not None and ub is not None:
        if lb == ub:
            return "E", lb
        raise ValueError(f"Constraint {name} is a range [{lb}, {ub}]; MATLAB needs one sense per row.")
    if ub is not None:
        return "L", ub
    if lb is not None:
        return "G", lb
    raise ValueError(f"Constraint {name} has no finite bound.")


def _variable_index(reactions: Sequence) -> Dict[str, Tuple[int, int]]:
    """Map solver variable name -> (reaction column, sign). Net flux = forward - reverse."""
    index = {}
    for j, rxn in enumerate(reactions):
        index[rxn.forward_variable.name] = (j, 1)
        index[rxn.reverse_variable.name] = (j, -1)
    return index


def _coupling_row(constraint, var_index, name: str) -> Dict[int, float]:
    """Read one coupling constraint back as {reaction column: coefficient on net flux}.

    Each reaction's forward and reverse variables must carry opposite coefficients;
    otherwise the constraint is not a constraint on net flux and has no C-row form.
    """
    coeffs = constraint.get_linear_coefficients(constraint.variables)
    by_rxn: Dict[int, Dict[int, float]] = {}
    for var, coef in coeffs.items():
        if coef == 0:
            continue
        if var.name not in var_index:
            raise ValueError(f"Constraint {name} uses variable {var.name}, which is not a reaction flux.")
        j, sign = var_index[var.name]
        by_rxn.setdefault(j, {})[sign] = float(coef)
    row = {}
    for j, parts in by_rxn.items():
        fwd, rev = parts.get(1, 0.0), parts.get(-1, 0.0)
        if not np.isclose(fwd, -rev):
            raise ValueError(
                f"Constraint {name}: forward/reverse coefficients {fwd}, {rev} for column {j} "
                "are not opposite, so it is not a constraint on net flux.")
        if fwd != 0:
            row[j] = fwd
    return row


def build_opt_problem_from_model(model) -> LPProblem:
    """Return the stacked LP [S; C] that MATLAB's buildOptProblemFromModel would build.

    MATLAB: LPproblem = buildOptProblemFromModel(model, verify, param)
    """
    notes = getattr(model, "notes", {}) or {}
    if notes.get("evars"):
        raise NotImplementedError("Extra variables (E/D) are not supported.")

    reactions = list(model.reactions)
    metabolites = list(model.metabolites)
    rxn_ids = [r.id for r in reactions]
    met_ids = [m.id for m in metabolites]
    met_row = {m.id: i for i, m in enumerate(metabolites)}
    n, m = len(reactions), len(metabolites)

    # Objective, bounds and S, column by column.
    c = np.zeros(n)
    lb = np.empty(n)
    ub = np.empty(n)
    rows, cols, vals = [], [], []
    for j, rxn in enumerate(reactions):
        c[j] = rxn.objective_coefficient
        lb[j], ub[j] = rxn.lower_bound, rxn.upper_bound
        for met, coef in rxn.metabolites.items():
            rows.append(met_row[met.id])
            cols.append(j)
            vals.append(coef)

    # Metabolite rows: sense and right-hand side from each mass-balance constraint
    # (cobrapy names it after the metabolite id). Usually all 'E' with b = 0.
    b = np.zeros(m)
    csense = np.empty(m, dtype="<U1")
    for i, met in enumerate(metabolites):
        con = model.constraints[met.id]
        csense[i], b[i] = _row_sense(con.lb, con.ub, met.id)

    # Coupling rows, in the order recorded when they were added.
    ctr_ids = list(notes.get(COUPLING_KEY, []))
    var_index = _variable_index(reactions) if ctr_ids else {}
    d = np.zeros(len(ctr_ids))
    dsense = np.empty(len(ctr_ids), dtype="<U1")
    for k, name in enumerate(ctr_ids):
        con = model.constraints[name]
        dsense[k], d[k] = _row_sense(con.lb, con.ub, name)
        for j, coef in _coupling_row(con, var_index, name).items():
            rows.append(m + k)
            cols.append(j)
            vals.append(coef)

    A = sp.csr_matrix((vals, (rows, cols)), shape=(m + len(ctr_ids), n))
    osense = -1 if model.objective.direction == "max" else 1

    return LPProblem(A=A, b=np.concatenate([b, d]), c=c, lb=lb, ub=ub,
                     csense=np.concatenate([csense, dsense]), osense=osense,
                     rxn_ids=rxn_ids, met_ids=met_ids, ctr_ids=ctr_ids)


def compare_to_matlab(lp: LPProblem, mat: dict, atol: float = 1e-9) -> List[str]:
    """Differences between a Python LPProblem and a MATLAB LPproblem loaded with scipy.io.loadmat.

    Export from MATLAB with:  LPproblem = buildOptProblemFromModel(model); save('lp.mat','LPproblem')
    and load with:            scipy.io.loadmat('lp.mat', squeeze_me=True)['LPproblem']
    Rows and columns are assumed to be in the same order (both come from the same .mat model).
    Returns an empty list when the two problems match.
    """
    def get(name):
        return mat[name] if isinstance(mat, dict) else getattr(mat, name)

    problems = []
    A_mat = sp.csr_matrix(get("A"))
    if A_mat.shape != lp.A.shape:
        return [f"A shape {lp.A.shape} vs MATLAB {A_mat.shape}"]
    diff = abs(lp.A - A_mat)
    if diff.nnz and diff.max() > atol:
        problems.append(f"A differs in {(diff > atol).nnz} entries")
    for name in ("b", "c", "lb", "ub"):
        py, ml = getattr(lp, name), np.ravel(get(name)).astype(float)
        bad = ~np.isclose(py, ml, atol=atol, equal_nan=True)
        if bad.any():
            problems.append(f"{name} differs at {int(bad.sum())} positions (first: {int(np.argmax(bad))})")
    ml_sense = np.array(list("".join(np.ravel(get("csense")).astype(str))))
    if ml_sense.shape != lp.csense.shape or (ml_sense != lp.csense).any():
        problems.append("csense differs")
    if int(get("osense")) != lp.osense:
        problems.append(f"osense {lp.osense} vs MATLAB {int(get('osense'))}")
    return problems
