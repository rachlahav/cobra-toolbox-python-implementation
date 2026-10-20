# Exported from the PSCM workbench (record 'changeCobraSolver', status: Drafting, updated 2026-10-05T21:00:00Z).
# Edit in the workbench, then re-run tools/export_workbench.py.
"""Select the solver used for each problem type (LP, MILP, QP, MIQP, ...).

MATLAB sources:
  src/base/solvers/changeCobraSolver.m
  src/base/solvers/getSetSolver/getCobraSolver.m
Translation status: Draft

`
------------------------------------------------------
solve_cobra_lp / qp / milp / miqp take an optProblem struct and solve it with
the solver chosen for that problem type,
obrapy's Model.solver only sets one solver for one mode (Gurobi)

cobrapy / optlang used
----------------------
cobra.util.solver.check_solver(name) returns the optlang interface module or
raises if that solver is not installed (imported by cobra.core.model). Solving
itself happens in the solve_cobra_* translations.
------------------------------------
Deviations from MATLAB
----------------------
1. Solver installation is checked by loading the optlang interface instead of
   the local MATLAB paths. 
2. TOMLAB and IBM CPLEX variants all map to 'cplex'. MATLAB-only solvers
   (pdco, minos, lindo, lp_solve, qpng, gurobi_mex, opti, matlab) are listed as
   not installed. 
3. validation_level > 0 calls solve_cobra_<type> from pscm.core.solve, which is
   not translated yet. Until it is, use validation_level 0 or -1.

-----------------------------
- validation_level == -1 sets the global before problem_type is upper-cased,
  so 'lp' sets nothing (solver_ok is NaN either way).
- Problem-type checks use contains(), a substring test: 'MIQP' passes a check
  against a solver whose types include 'QP'.
- The test problem omits F when solver is 'mosek' with type 'CLP' or when type
  is 'all'; 'all' never reaches that line.
"""
import math
from typing import Dict, Optional, Tuple

from cobra.util.solver import check_solver

# MATLAB: global OPT_PROB_TYPES
OPT_PROB_TYPES = ["LP", "MILP", "QP", "MIQP", "NLP", "EP", "CLP"]

# MATLAB: global CBT_<type>_SOLVER
CBT_SOLVER: Dict[str, Optional[str]] = {t: None for t in OPT_PROB_TYPES}

# MATLAB solver name -> optlang interface name (None = not available in Python)
SOLVER_INTERFACE = {
    "gurobi": "gurobi",
    "ibm_cplex": "cplex", "cplexlp": "cplex", "tomlab_cplex": "cplex", "cplex_direct": "cplex",
    "glpk": "glpk",
    "mosek": None, "pdco": None, "quadMinos": None, "dqqMinos": None, "matlab": None,
    "lindo_new": None, "lindo_legacy": None, "lp_solve": None, "qpng": None,
    "gurobi_mex": None, "opti": None, "tomlab_snopt": None,
    # Python-only, through optlang
    "scipy": "scipy", "osqp": "osqp", "hybrid": "hybrid", "coinor_cbc": "coinor_cbc",
}

# Problem types each solver can solve (from the changeCobraSolver.m header;
# Python-only solvers from their optlang capabilities).
SOLVER_TYPES = {
    "gurobi": ["LP", "MILP", "QP", "MIQP"],
    "ibm_cplex": ["LP", "MILP", "QP", "MIQP"], "tomlab_cplex": ["LP", "MILP", "QP", "MIQP"],
    "cplex_direct": ["LP", "MILP", "QP", "MIQP"], "cplexlp": ["LP"],
    "glpk": ["LP", "MILP"],
    "mosek": ["LP", "MILP", "QP", "EP", "CLP"], "pdco": ["LP", "QP"],
    "quadMinos": ["LP", "NLP"], "dqqMinos": ["LP", "QP"], "matlab": ["LP", "NLP"],
    "lindo_new": ["LP"], "lindo_legacy": ["LP"], "lp_solve": ["LP"], "qpng": ["QP"],
    "gurobi_mex": ["LP", "MILP", "QP", "MIQP"], "opti": ["LP", "MILP", "QP"],
    "tomlab_snopt": ["NLP"],
    "scipy": ["LP"], "osqp": ["LP", "QP"], "hybrid": ["LP", "QP"], "coinor_cbc": ["LP", "MILP"],
}

# MATLAB: global SOLVERS (fields installed, working, type). Built on first use.
SOLVERS: Dict[str, dict] = {}


def _build_solvers() -> None:
    """MATLAB initCobraToolbox fills SOLVERS; here: try to load each optlang interface."""
    for name, iface in SOLVER_INTERFACE.items():
        installed = False
        if iface is not None:
            try:
                check_solver(iface)
                installed = True
            except Exception:
                installed = False
        SOLVERS[name] = {"installed": installed, "working": installed,
                         "type": SOLVER_TYPES[name], "interface": iface}


def _contains(text: str, patterns) -> bool:
    """MATLAB contains(str, cellOfPatterns): true if any pattern is a substring."""
    return any(p in text for p in patterns)


def change_cobra_solver(solver_name: Optional[str] = None, problem_type: Optional[str] = None,
                        print_level: int = 1, validation_level: int = 1) -> Tuple[object, bool]:
    """Change the solver used for a problem type.

    MATLAB: [solverOK, solverInstalled] = changeCobraSolver(solverName, problemType,
                                                            printLevel, validationLevel)

    validation_level: -1 assign only; 0 check installed, don't validate;
                      1 validate silently (default); 2 validate and show output.
    Returns (solver_ok, solver_installed); solver_ok is NaN when not validated (-1).
    """
    if solver_name == "cplex":
        solver_name = "ibm_cplex"

    solver_installed = True

    # --- assign only (MATLAB: validationLevel == -1) ---
    if validation_level == -1:
        if problem_type in CBT_SOLVER:          # quirk kept: not upper-cased yet
            CBT_SOLVER[problem_type] = solver_name
        return math.nan, solver_installed

    solver_installed = False

    if not SOLVERS:                              # MATLAB: isempty(SOLVERS) -> initCobraToolbox
        _build_solvers()

    # --- no arguments: print defined solvers (MATLAB: nargin < 1) ---
    if solver_name is None:
        defined = {t: s for t, s in CBT_SOLVER.items() if s}
        if not defined:
            print("No solvers are defined!")
        else:
            print("Defined solvers are:")
            for t in OPT_PROB_TYPES:
                if CBT_SOLVER[t]:
                    print(f"    CBT_{t}_SOLVER: {CBT_SOLVER[t]}")
        return None, solver_installed

    # legacy support for other versions of gurobi
    if solver_name.lower() in ("gurobi6", "gurobi7"):
        solver_name = "gurobi"

    # legacy support for MPS
    if solver_name.lower() == "mps":
        print(" > changeCobraSolver: The interface to 'mps' from 'changeCobraSolver()' is no longer supported.")
        raise ValueError(" -> Use writeCbModel(model, 'mps') instead.")

    problem_type = "LP" if problem_type is None else problem_type.upper()

    # --- solver must be known ---
    if solver_name not in SOLVERS:
        raise ValueError(f"changeCobraSolver: The solver {solver_name} is not available to the "
                         "COBRA Toolbox on your machine.")
    if validation_level == 0:
        if not SOLVERS[solver_name]["installed"]:
            if print_level > 0:
                print(f" > changeCobraSolver: Solver {solver_name} is not installed.")
            return False, SOLVERS[solver_name]["installed"]
        if not SOLVERS[solver_name]["working"]:
            if print_level > 0:
                print(f" > changeCobraSolver: Solver {solver_name} is installed but not working properly.")
            return False, SOLVERS[solver_name]["installed"]

    # --- 'all': set every problem type this solver supports ---
    if problem_type == "ALL":                    # MATLAB: strcmpi(problemType, 'all')
        solved_problems = SOLVERS[solver_name]["type"]
        solver_ok = False
        for t in solved_problems:
            solver_ok, solver_installed = change_cobra_solver(solver_name, t, print_level, validation_level)
            if print_level > 0:
                print(f" > changeCobraSolver: Solver for {t} problems has been set to {solver_name}.")
        for t in [t for t in OPT_PROB_TYPES if t not in solved_problems]:
            used = CBT_SOLVER[t]
            info = "No solver set for this problemtype" if not used else f"Currently used: {used}"
            if print_level > 0:
                print(f" > changeCobraSolver: Solver {solver_name} not supported for problems of type {t}. {info}")
        return solver_ok, solver_installed

    # --- problem type checks (MATLAB contains(): substring test, quirk kept) ---
    solver_ok = False
    if not _contains(problem_type, OPT_PROB_TYPES):
        solver_installed = SOLVERS[solver_name]["installed"]
        if print_level > 0:
            raise ValueError(f"changeCobraSolver: {problem_type}  problems cannot be solved in The COBRA Toolbox")
        return solver_ok, solver_installed
    if not _contains(problem_type, SOLVERS[solver_name]["type"]):
        solver_installed = SOLVERS[solver_name]["installed"]
        if print_level > 0:
            raise ValueError(f"Solver {solver_name} cannot solve {problem_type} problems")
        return solver_ok, solver_installed

    # --- installation check (MATLAB: checkSolverInstallationFile etc.) ---
    if SOLVERS[solver_name]["interface"] is None:
        if print_level > 0:
            raise ValueError(f"changeCobraSolver: Solver {solver_name} is not available in Python "
                             "(no optlang interface).")
        solver_ok = False
    else:
        solver_ok = SOLVERS[solver_name]["installed"]
        if not solver_ok and print_level > 0:
            raise ValueError(f"changeCobraSolver: Solver {solver_name} is not installed.")

    # --- set the global, validating with a tiny problem ---
    if solver_ok:
        solver_installed = True
        if validation_level > 0:
            old_val = CBT_SOLVER.get(problem_type)
            CBT_SOLVER[problem_type] = solver_name
            # MATLAB test problem: A=[0 1], b=0, c=[1;1], osense=-1, lb=ub=0, csense='E',
            # vartype ['C';'I'], x0=[0;0], F=speye(2) (omitted for mosek CLP / 'all')
            problem = {"A": [[0, 1]], "b": [0], "c": [1, 1], "osense": -1,
                       "lb": [0, 0], "ub": [0, 0], "csense": "E", "vartype": ["C", "I"],
                       "x0": [0, 0], "names": []}
            if not ((solver_name == "mosek" and problem_type == "CLP") or problem_type == "all"):
                problem["F"] = [[1, 0], [0, 1]]
            try:
                if problem_type != "CLP":       # skip the CLP solver until further developments
                    from pscm.core import solve  # solve_cobra_lp / qp / milp / miqp (Option B)
                    solve_fn = getattr(solve, f"solve_cobra_{problem_type.lower()}")
                    solve_fn(problem, printLevel=3 if validation_level > 1 else 0)
            except Exception as err:
                if print_level > 0:
                    print(f"There was an error! The message was:\n{err}")
                    raise
                solver_ok = False
                CBT_SOLVER[problem_type] = old_val
        else:
            # if unvalidated, simply set the solver without testing
            CBT_SOLVER[problem_type] = solver_name
    else:
        if solver_name == "gurobi":
            print("Gurobi installed, and licence current? Check GUROBI_HOME and GRB_LICENSE_FILE.")
        else:
            print(f"Could not find installation of {solver_name}, so it cannot be tested")

    return solver_ok, solver_installed


def get_cobra_solver(solver_type: str = "LP", validate: bool = False) -> Tuple[Optional[str], object]:
    """Return the solver currently set for a problem type.

    MATLAB: [solverName, solverOK] = getCobraSolver(solverType, validate)
    solver_ok is NaN unless validate is true.
    """
    solver_ok = math.nan
    solver_name = CBT_SOLVER.get(solver_type)    # MATLAB: switch on 'LP','QP','MILP','MIQP','NLP'
    if validate:
        solver_ok, _ = change_cobra_solver(solver_name, solver_type, 0, 1)
        if not solver_ok:
            # if there is something wrong, do a verbose validation
            solver_ok, _ = change_cobra_solver(solver_name, solver_type, 1, 1)
    return solver_name, solver_ok
