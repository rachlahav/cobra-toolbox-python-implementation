# Exported from the PSCM workbench (record 'optimizeCbModel', status: Drafting, updated 2026-10-12T23:26:17.819Z).
# Edit in the workbench, then re-run tools/export_workbench.py.
"""Flux balance analysis with optional minimisation of the one or two norm, on cobrapy's core.
    solveCobraLP / solveCobraQP          -> model.optimize()
    row c'v = objectiveLP  ([A; c'])     -> a constraint on the linear objective
    1-norm with p, q >= 0 (v = p - q)    -> cobrapy's forward/reverse variables
    F = diag(minNorm)                    -> quadratic objective on reaction.flux_expression
    addLoopLawConstraints                -> not available (MATLAB source not in project)
    whole-body 'QRLP' problem            -> extra optlang variables r, p, z

"""
import time
import warnings

import numpy as np
from optlang.symbolics import Zero, add


from pscm.core.solver_params import get_cobra_solver_params

_NOT_GIVEN = object()                    # MATLAB: argument not passed

# optlang interfaces that can solve QPs / MIQPs (as listed in changeCobraSolver.m
# for the solvers that have an optlang interface)
_QP_SOLVERS = {"gurobi", "cplex", "osqp", "hybrid"}
_MIQP_SOLVERS = {"gurobi", "cplex"}
COUPLING_KEY = "coupling_constraint_ids"
MATLAB_FIELDS_KEY = "matlab_fields"


def optimize_cb_model(model, osense_str=_NOT_GIVEN, min_norm=_NOT_GIVEN,
                      allow_loops=_NOT_GIVEN, param=_NOT_GIVEN) -> dict:
    """MATLAB: solution = optimizeCbModel(model, osenseStr, minNorm, allowLoops, param)

    model       cobrapy Model (whole-body conventions: coupling constraints and
                MATLAB-only fields in model.notes)
    osense_str  'max' / 'min', or a param dict (then minNorm/allowLoops come from it)
    min_norm    None/0 (FBA only), 'one', 'oneInternal', a positive scalar or an
                nRxns vector (2-norm), True (= 1e-6)
    allow_loops True (default); False (loopless) raises until addLoopLawConstraints is translated
    param       dict: minNorm, printLevel, feasTol, optTol, timeLimit,
                solveWBMmethod ('LP', 'QP', 'QRLP'), solver (optlang name), solver-native settings

    Returns the MATLAB solution struct as a dict.
    """
    fields = dict((getattr(model, "notes", {}) or {}).get(MATLAB_FIELDS_KEY, {}))
    model_osense_str = model.objective.direction             # MATLAB model.osenseStr

    if osense_str is not _NOT_GIVEN:
        if _isempty(osense_str):
            model_osense_str = "max"
        else:
            if isinstance(osense_str, dict):
                if param is not _NOT_GIVEN:
                    raise ValueError("osenseStr is a structure and param structure is also present")
                else:
                    param = osense_str
                    # quirk kept: a cobrapy model always "has" osenseStr, so this is always 'max'
                    model_osense_str = "max"
                    min_norm = param["minNorm"] if "minNorm" in param else None
                    allow_loops = param["allowLoops"] if "allowLoops" in param else 1
            else:
                model_osense_str = osense_str

    if param is _NOT_GIVEN:
        param = {}

    # --- lines 266-294: minNorm ---
    if min_norm is not _NOT_GIVEN:
        if isinstance(min_norm, bool):                     # islogical
            min_norm = 1e-6 if min_norm else None
        if _isnumeric(min_norm) and np.ndim(min_norm) == 0 and min_norm == 0:   # isequal(minNorm,0)
            min_norm = None
        if _strcmp(min_norm, "zero"):
            raise NotImplementedError("minNorm = 'zero' needs sparseLP, which is not translated yet.")
    else:
        min_norm, = get_cobra_solver_params("LP", "minNorm")

    # --- lines 296-302: allowLoops ---
    if allow_loops is _NOT_GIVEN or _isempty(allow_loops):
        allow_loops = True

    # --- line 304: solver parameters ---
    print_level, primal_only_flag, verify, feas_tol, _ = get_cobra_solver_params(
        "LP", ["printLevel", "primalOnly", "verify", "feasTol"], param)

    # --- lines 306-322: dimensions ---
    reactions = list(model.reactions)
    metabolites = list(model.metabolites)
    n_mets, n_rxns = len(metabolites), len(reactions)
    ctr_ids = list((getattr(model, "notes", {}) or {}).get(COUPLING_KEY, []))
    model_c = 1 if ctr_ids else 0
    n_ctrs = len(ctr_ids)
    if "E" in fields:
        raise NotImplementedError("Extra variables (model.E) are not supported.")
    model_e, n_vars = 0, 0
    row_ids = [m.id for m in metabolites] + ctr_ids          # rows of [S; C]

    s_consistent_rxn_bool = None
    if _strcmp(min_norm, "oneInternal"):
        s_consistent_rxn_bool = np.asarray(fields["SConsistentRxnBool"], dtype=bool)

    # --- lines 330-338: legacy model.g ---
    if "g" in fields and "cf" not in fields:
        if "g0" in fields or "g1" in fields:
            warnings.warn("model.g ignored by optimizeCbModel. zero and one norm weights are "
                          "separately specified in model.g0 and model.g1 respectively")
        else:
            raise ValueError("model.g no longer supported by optimizeCbModel. zero and one norm weights "
                             "must be separately specified in model.g0 and model.g1 respectively")

    # --- lines 340-380: norm weights ---
    zero_norm_weights = None
    if "g0" in fields:
        if len(fields["g0"]) != n_rxns:
            raise ValueError("model.g0 must be nRxns x 1")
        zero_norm_weights = np.ravel(fields["g0"]).astype(float)
    one_norm_weights = None
    if "g1" in fields:
        if len(fields["g1"]) != n_rxns:
            raise ValueError("model.g1 must be nRxns x 1")
        one_norm_weights = np.ravel(fields["g1"]).astype(float)
    if not _isempty(zero_norm_weights):
        if np.all(zero_norm_weights == 1) and isinstance(min_norm, str):
            if _strcmp("minNorm", "optimizeCardinality") and not _isempty(one_norm_weights):
                min_norm = "one"                            # quirk kept: never reached
    if "g2" in fields:
        raise ValueError("minimisation of two norm in combination with one and zero norm is not yet supported")

    # MATLAB model.c, model.lb, model.ub
    c = np.array([r.objective_coefficient for r in reactions], dtype=float)
    lb = np.array([r.lower_bound for r in reactions], dtype=float)
    ub = np.array([r.upper_bound for r in reactions], dtype=float)

    # --- lines 384-398: decide whether to solve the LP first ---
    do_linear_optimisation_first = 1
    if np.all(c == 0):
        do_linear_optimisation_first = 0
    if _isempty(min_norm):
        do_linear_optimisation_first = 1
    if np.all(((lb == ub) & (c != 0)) == (c != 0)) and not _isempty(min_norm):
        do_linear_optimisation_first = 0

    wbm_method = param.get("solveWBMmethod")
    with model:                                              # every change below is undone on exit
        model.objective.direction = model_osense_str
        osense = -1 if model_osense_str == "max" else 1

        # --- lines 400-407: whole-body QP methods ---
        if wbm_method in ("QP", "QRLP", "QRQP"):
            c = np.zeros(n_rxns)                             # model.c(:)=0
            _set_objective(model, direction=model_osense_str)
            do_linear_optimisation_first = 0
            min_norm = param["minNorm"]
        if wbm_method == "QRQP":
            raise NotImplementedError("solveWBMmethod 'QRQP': its buildOptProblemFromModel section is not yet read.")

        # Solver: MATLAB param.solver, else the global CBT_*_SOLVER. cobrapy: the model's own
        # solver; param["solver"] is set on the model (reverted on exit). Done before any
        # variable is added, because changing model.solver clones the problem.
        if "solver" in param and param["solver"]:
            model.solver = param["solver"]
        needs_qp = (_isnumeric(min_norm) and (np.size(min_norm) > 1 or np.all(np.asarray(min_norm) > 0))) \
            or _strcmp(min_norm, "1e-6")
        if needs_qp and _solver_name(model) not in _QP_SOLVERS:
            raise ValueError(f"The model's solver '{_solver_name(model)}' cannot solve QPs (minNorm > 0). "
                             "Set model.solver = 'gurobi' or 'cplex', or pass param['solver'].")

        # --- line 409 (buildOptProblemFromModel) -> whole-body QRLP reformulation ---
        extra_vars = []
        if wbm_method == "QRLP":
            # MATLAB: variables [v r p z]; A v + r = b; v - p - z = 0; v, r, p free; lb <= z <= ub
            r_vars = [model.problem.Variable(f"QRLP_r_{i}") for i in range(len(row_ids))]
            p_vars = [model.problem.Variable(f"QRLP_p_{j}") for j in range(n_rxns)]
            z_vars = [model.problem.Variable(f"QRLP_z_{j}", lb=lb[j], ub=ub[j]) for j in range(n_rxns)]
            model.add_cons_vars(r_vars + p_vars + z_vars)
            model.solver.update()
            for i, row in enumerate(row_ids):
                model.constraints[row].set_linear_coefficients({r_vars[i]: 1.0})
            link = []
            for j, rxn in enumerate(reactions):
                rxn.bounds = (-np.inf, np.inf)               # v is free; its bounds move to z
                con = model.problem.Constraint(Zero, lb=0, ub=0, name=f"QRLP_link_{j}")
                link.append(con)
            model.add_cons_vars(link)
            model.solver.update()
            for j, rxn in enumerate(reactions):
                link[j].set_linear_coefficients({rxn.forward_variable: 1.0, rxn.reverse_variable: -1.0,
                                                 p_vars[j]: -1.0, z_vars[j]: -1.0})
            extra_vars = [v.name for v in r_vars + p_vars + z_vars]     # looked up by name later
            # (MATLAB's F = 1e-16 / max|c|*100 regularisation is overwritten by F = diag(minNorm) below)
        n_total_vars = n_rxns + len(extra_vars)

        t1 = time.perf_counter()
        solution = None
        objective_lp = None

        # --- lines 418-497: linear optimisation first ---
        if do_linear_optimisation_first:
            if allow_loops:
                param_lp = dict(param)
                param_lp.pop("minNorm", None)
                solution = _solve(model, reactions, row_ids, osense, c, param_lp, "LP")
            else:
                # MATLAB: MILPproblem = addLoopLawConstraints(optProblem, model, 1:nRxns);
                #         solution = solveCobraMILP(MILPproblem);
                _add_loop_law_constraints(model, np.arange(n_rxns))
                solution = _solve(model, reactions, row_ids, osense, c, param, "MILP")
            objective_lp = solution["obj"]

        f_on_opt_problem = False      # MATLAB: is F stored on optProblem (not optProblem2)?
        F_diag = None

        # --- lines 499-866: second step, minimise a norm ---
        if (do_linear_optimisation_first == 0 and not _isempty(min_norm)) or \
           (do_linear_optimisation_first == 1 and solution["stat"] == 1 and not _isempty(min_norm)):

            if _strcmp(min_norm, "optimizeCardinality"):
                raise NotImplementedError("minNorm = 'optimizeCardinality' needs optimizeCardinality, "
                                          "which is not translated yet.")

            elif _strcmp(min_norm, "one") or _strcmp(min_norm, "oneInternal"):
                # 'one': min sum w|v| over all reactions; 'oneInternal': only SConsistentRxnBool.
                # MATLAB p, q >= 0 are cobrapy's forward / reverse variables.
                if _strcmp(min_norm, "one"):
                    selected = np.arange(n_rxns)
                    weights = one_norm_weights if not _isempty(one_norm_weights) else np.ones(n_rxns)
                    if _isempty(objective_lp):
                        objective_lp = 0
                else:
                    selected = np.flatnonzero(s_consistent_rxn_bool)
                    weights = np.ones(n_rxns)
                    if not _isempty(one_norm_weights):
                        weights = one_norm_weights
                        if np.any(one_norm_weights[selected] < 0):
                            warnings.warn("minNorm = 'oneInternal' may not eliminate thermodynamically infeasible "
                                          "fluxes if model.g1(SConsistentRxnBool) entries are negative")
                    if _isempty(objective_lp):
                        raise ValueError("oneInternal: objectiveLP is empty (MATLAB fails on matrix dimensions)")
                # objective row: c'v >= objectiveLP when maximising, <= when minimising
                if osense == -1:
                    _linear_objective_constraint(model, reactions, c, lb=objective_lp)
                else:
                    _linear_objective_constraint(model, reactions, c, ub=objective_lp)
                linear = {}
                for j in selected:
                    linear[reactions[j].forward_variable] = float(weights[j])
                    linear[reactions[j].reverse_variable] = float(weights[j])
                _set_objective(model, linear=linear, direction="min")      # optProblem2.osense = 1
                if not allow_loops and _strcmp(min_norm, "one"):
                    # MATLAB: addLoopLawConstraints(optProblem2, model, 1:nRxns) + solveCobraMILP
                    _add_loop_law_constraints(model, np.arange(n_rxns))
                solution = _solve(model, reactions, row_ids, osense, np.zeros(n_rxns), param,
                                  "LP" if allow_loops else "MILP")

            elif _strcmp(min_norm, "QRLP"):
                raise NotImplementedError("minNorm = 'QRLP' (as a string) is not supported; use "
                                          "param.solveWBMmethod = 'QRLP' with a numeric minNorm.")

            elif _strcmp(min_norm, "QRQP"):
                # MATLAB: minNorm = minNormWBM;  -- undefined, so MATLAB errors here
                raise NameError("Undefined function or variable 'minNormWBM' (as in optimizeCbModel.m)")

            elif (isinstance(min_norm, str) and len(min_norm) >= 1) or \
                 (_isnumeric(min_norm) and (np.size(min_norm) > 1 or np.all(np.asarray(min_norm) > 0))):
                # --- 2-norm: min 0.5 * sum(minNorm_i * x_i^2) ---
                if _isnumeric(min_norm):
                    arr = np.asarray(min_norm, dtype=float)
                    if arr.size == n_total_vars and arr.ndim == 1:
                        min_norm = arr.ravel()
                    elif arr.size == 1:
                        min_norm = np.ones(n_total_vars) * float(arr)
                    else:
                        raise ValueError(f"minNorm has dimensions {arr.shape} but it can only of the form "
                                         "{(0), 'one', 'zero', > 0 , n x 1 vector}.")
                elif min_norm == "1e-6":
                    min_norm = np.ones(n_total_vars) * 1e-6
                else:
                    raise ValueError(f"minNorm has dimensions 1 x {len(min_norm)} but it can only of the form "
                                     "{(0), 'one', 'zero', > 0 , n x 1 vector}.")
                F_diag = min_norm
                x_exprs = [r.flux_expression for r in reactions] + [model.variables[n] for n in extra_vars]
                quad = [0.5 * F_diag[i] * x_exprs[i] ** 2 for i in range(n_total_vars) if F_diag[i] != 0]

                if do_linear_optimisation_first:
                    # optProblem2: add row c'v = objectiveLP, osense = 1, objective c'x + 0.5 x'Fx
                    _linear_objective_constraint(model, reactions, c, lb=objective_lp, ub=objective_lp)
                    lin_sign = 1
                else:
                    # optProblem keeps its osense: objective osense*c'x + 0.5 x'Fx
                    f_on_opt_problem = True
                    lin_sign = osense
                linear = {}
                for j in np.flatnonzero(c):
                    linear[reactions[j].forward_variable] = lin_sign * float(c[j])
                    linear[reactions[j].reverse_variable] = -lin_sign * float(c[j])
                _set_objective(model, linear=linear, quadratic=quad, direction="min")   # solveCobraQP minimises
                if not allow_loops:
                    # MATLAB: addLoopLawConstraints(optProblem(2), model, 1:nTotalVars) + solveCobraMIQP
                    # (quirk: 1:nTotalVars, not 1:nRxns)
                    if _solver_name(model) not in _MIQP_SOLVERS:
                        raise ValueError(f"The model's solver '{_solver_name(model)}' cannot solve MIQPs.")
                    _add_loop_law_constraints(model, np.arange(n_total_vars))
                solution = _solve(model, reactions, row_ids, osense, c, param,
                                  "QP" if allow_loops else "MIQP", extra_vars=extra_vars, F_diag=F_diag)

    if solution is None:
        raise ValueError("No problem was solved (no linear step and no valid minNorm).")

    # --- lines 868-897: report status ---
    solution["f0"] = np.nan
    solution["f1"] = np.nan
    solution["f2"] = np.nan
    stat = solution["stat"]
    if stat == 1:
        if print_level > 0:
            print("Optimal solution found.")
    elif stat == -1:
        if print_level > 0:
            warnings.warn("No solution reported (timelimit, numerical problem etc).")
    elif stat == 0:
        if print_level > 0:
            warnings.warn("Infeasible model.")
    elif stat == 2:
        if print_level > 0:
            warnings.warn("Unbounded model.")
    elif stat == 3:
        if print_level > 0:
            warnings.warn("Solution exists, but either scaling problems or not proven to be optimal.")
    else:
        print(stat)
        raise ValueError("solution.stat must be in {-1, 0 , 1, 2, 3}")

    if "dual" not in solution or _isempty(solution["dual"]):
        primal_only_flag = 1

    # --- lines 905-1027: unpack an optimal solution ---
    if solution["stat"] == 1 or solution["stat"] == 3:
        full = solution["full"]
        if _isempty(objective_lp):
            solution["f"] = objective_lp                   # quirk kept: None when no linear step
        else:
            solution["f"] = float(np.dot(c, full[:n_rxns]))

        if _isempty(min_norm):
            min_norm = "empty"
        if _isnumeric(min_norm):
            min_norm = "two"

        if min_norm == "empty":
            solution["f1"] = float(np.dot(c, full[:n_rxns]))
        elif min_norm == "one":
            solution["f1"] = float(np.sum(np.abs(full[:n_total_vars])))
        elif min_norm == "two":
            solution["f1"] = float(np.dot(c, full[:n_rxns]))
            if f_on_opt_problem:
                solution["f2"] = 0.5 * float(np.dot(F_diag, full[:n_total_vars] ** 2))
            else:
                print(param)
                warnings.warn("optProblem.F missing")    # quirk kept: F was on optProblem2
        # other norms ('oneInternal'): quirk kept, f0 = f1 = f2 = NaN

        solution["obj"] = float(np.nansum([solution["f0"], solution["f1"], solution["f2"]]))
        solution["v"] = full[:n_rxns]
        solution["vars_v"] = None                         # model.E not supported
        solution["x"] = solution["v"]

        if not _isempty(solution.get("dual")):
            solution["y"] = solution["dual"][:n_mets]
            if model_c:
                solution["ctrs_y"] = solution["dual"][n_mets:n_mets + n_ctrs]
        if not _isempty(solution.get("rcost")):
            solution["w"] = solution["rcost"][:n_rxns]
        if not _isempty(solution.get("slack")):
            solution["s"] = solution["slack"][:n_mets]
            if model_c:
                solution["ctrs_s"] = solution["slack"][n_mets:n_mets + n_ctrs]

        if not primal_only_flag and allow_loops:
            solution["y"] = solution["dual"][:n_mets]
            if model_c:
                solution["ctrs_y"] = solution["dual"][n_mets:n_mets + n_ctrs]
            solution["w"] = solution["rcost"]             # quirk kept: full vector again
            solution["s"] = solution["slack"]             # quirk kept: full vector again

        solution["time"] = time.perf_counter() - t1

        field_order = ["f", "f0", "f1", "f2", "v", "y", "w", "s", "solver", "lpmethod", "qpmethod",
                       "stat", "origStat", "time", "basis", "vars_v", "vars_w", "ctrs_y", "ctrs_s",
                       "x", "full", "obj", "rcost", "dual", "slack"]
        absent = [k for k in solution if k not in field_order]
        present = [k for k in field_order if k in solution]
        solution = {k: solution[k] for k in absent + present}
    else:
        solution.update(obj=np.nan, f=np.nan, f0=np.nan, f1=np.nan, f2=np.nan,
                        v=None, y=None, w=None, s=None)
        if model_c:
            solution["ctrs_y"] = None
            solution["ctrs_s"] = None
        solution["x"] = solution["v"]
        solution["time"] = time.perf_counter() - t1

    # --- lines 1051-1059: remove fields not in the output specification ---
    for key in ("full", "dual", "rcost", "slack"):
        solution.pop(key, None)
    return solution



# ---------------------------------------------------------------- helpers ---
def _solver_name(model) -> str:
    """Name of the model's solver, e.g. 'gurobi' (from optlang's 'optlang.gurobi_interface')."""
    return model.problem.__name__.split(".")[-1].replace("_interface", "")


def _add_loop_law_constraints(*_):
    """MATLAB addLoopLawConstraints: not translated yet."""
    raise NotImplementedError(
        "allowLoops = 0 needs addLoopLawConstraints.m, which is not in the project files. "
        "Add its source (cobratoolbox src/analysis/...) to translate the loop-law MILP.")


def _isempty(x) -> bool:
    """MATLAB isempty for the values used here (None stands for [])."""
    if x is None:
        return True
    if isinstance(x, str):
        return len(x) == 0
    return np.size(x) == 0


def _strcmp(x, s: str) -> bool:
    """MATLAB strcmp: true only for a string equal to s."""
    return isinstance(x, str) and x == s


def _isnumeric(x) -> bool:
    if x is None or isinstance(x, (str, bool)):
        return False
    return np.issubdtype(np.asarray(x).dtype, np.number)


def _linear_objective_constraint(model, reactions, c, lb=None, ub=None, name="objectiveLP"):
    """MATLAB row c'v (csense) objectiveLP, as a cobrapy constraint on net fluxes."""
    con = model.problem.Constraint(Zero, lb=lb, ub=ub, name=name)
    model.add_cons_vars([con])
    model.solver.update()
    coeffs = {}
    for j in np.flatnonzero(c):
        coeffs[reactions[j].forward_variable] = float(c[j])
        coeffs[reactions[j].reverse_variable] = -float(c[j])
    con.set_linear_coefficients(coeffs)
    return con


def _set_objective(model, linear=None, quadratic=(), direction="min"):
    """Set objective  sum(linear) + sum(quadratic terms)  (quadratic as expressions)."""
    model.objective = model.problem.Objective(add(list(quadratic)) if quadratic else Zero,
                                              direction=direction, sloppy=True)
    if linear:
        model.solver.update()
        model.objective.set_linear_coefficients(linear)


def _apply_params(model, param, problem_type):
    """Apply COBRA solver parameters (and gurobi-native ones) to model.solver; return a restore function."""
    feas_tol, opt_tol, time_limit, print_level, _ = get_cobra_solver_params(
        problem_type, ["feasTol", "optTol", "timeLimit", "printLevel"], param)
    cfg = model.solver.configuration
    old = {"feas": cfg.tolerances.feasibility, "opt": cfg.tolerances.optimality,
           "timeout": cfg.timeout, "verbosity": cfg.verbosity}
    cfg.tolerances.feasibility = feas_tol
    cfg.tolerances.optimality = opt_tol
    cfg.timeout = None if time_limit >= 1e36 else time_limit
    cfg.verbosity = int(min(max(print_level, 0), 3))
    old_native = {}
    if _solver_name(model) == "gurobi":
        params = model.solver.problem.Params
        for key, value in param.items():
            if hasattr(params, key):
                old_native[key] = getattr(params, key)
                setattr(params, key, value)

    def restore():
        cfg.tolerances.feasibility = old["feas"]
        cfg.tolerances.optimality = old["opt"]
        cfg.timeout = old["timeout"]
        cfg.verbosity = old["verbosity"]
        if old_native:
            params = model.solver.problem.Params
            for key, value in old_native.items():
                setattr(params, key, value)
    return restore


def _solve(model, reactions, row_ids, osense, c, param, problem_type, extra_vars=(), F_diag=None):
    """Solve the model's current problem; return MATLAB's solveCobraLP/QP solution fields.

    full  = [net fluxes, extra variables]       rcost = osense * reduced costs (same order)
    dual  = osense * duals of [metabolite rows; coupling rows]
    slack = b - A*x for those rows             obj = c'x (+ 0.5 x'Fx for a QP)
    """
    restore = _apply_params(model, param, problem_type)
    t0 = time.perf_counter()
    try:
        sol = model.optimize(raise_error=False)          # cobra.core.model / get_solution
        orig_stat = sol.status
        if orig_stat == "optimal":
            stat = 1
        elif orig_stat == "infeasible":
            stat = 0
        elif orig_stat == "unbounded":
            stat = 2
        elif orig_stat == "infeasible_or_unbounded":
            # solveCobraLP: remove the objective and solve again
            with model:
                model.objective = model.problem.Objective(Zero, direction="min")
                stat = 2 if model.solver.optimize() == "optimal" else 0
        else:
            stat = -1
    finally:
        restore()
    out = {"solver": _solver_name(model), "stat": stat, "origStat": orig_stat,
           "time": time.perf_counter() - t0}
    if stat != 1:
        out.update(full=None, obj=np.nan, dual=None, rcost=None, slack=None)
        return out

    def _safe(get):
        try:
            value = get()
            return np.nan if value is None else value
        except Exception:                     # e.g. no duals for a MILP
            return np.nan

    extra = [model.variables[n] for n in extra_vars]          # extra_vars holds variable names
    full = np.concatenate([sol.fluxes.values, np.array([v.primal for v in extra])])
    rcost = np.concatenate([sol.reduced_costs.values,
                            np.array([_safe(lambda v=v: v.dual) for v in extra])])
    duals = np.array([_safe(lambda r=r: model.constraints[r].dual) for r in row_ids], dtype=float)
    slack = []
    for r in row_ids:
        con = model.constraints[r]
        rhs = con.lb if con.lb is not None else con.ub      # b for this row
        slack.append(rhs - con.primal)                       # A*x + s = b
    obj = float(np.dot(c, full[:len(c)]))
    if F_diag is not None:
        obj += 0.5 * float(np.dot(F_diag, full[:len(F_diag)] ** 2))
    out.update(
        full=full,
        obj=obj,
        dual=None if np.all(np.isnan(duals)) else osense * duals,
        rcost=None if np.all(np.isnan(rcost)) else osense * rcost,
        slack=np.array(slack),
    )
    return out
