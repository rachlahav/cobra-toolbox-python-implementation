# Exported from the PSCM workbench (record 'getCobraSolverParams', status: Drafting, updated 2026-10-13T00:05:47.186Z).
# Edit in the workbench, then re-run tools/export_workbench.py.
"""COBRA solver parameter defaults and lookup.

MATLAB sources
  src/base/solvers/param/getCobraSolverParams.m 
  src/base/solvers/param/getCobraSolverParamsOptionsForType.m

------------
Returns the value of each requested solver parameter, taking the first that is
set from: the param structure passed in > the global parameters for the
problem type > the defaults defined in the function. optimizeCbModel and
optimizeWBModel use it for minNorm, printLevel, primalOnly, verify and feasTol.

"""
from typing import Dict, List, Optional, Sequence, Union

PROBLEM_TYPES = ("LP", "QP", "MILP", "EP", "CLP", "MIQP", "NLP")

# MATLAB: global CBT_<type>_PARAMS
GLOBAL_SOLVER_PARAMS: Dict[str, dict] = {t: {} for t in PROBLEM_TYPES}

# MATLAB getCobraSolverParamsOptionsForType (LP and QP read from source)
_OPTIONS_FOR_TYPE = {
    "LP": ["multiscale", "problemType", "verify", "minNorm", "printLevel", "primalOnly",
           "saveInput", "feasTol", "optTol", "solver", "debug", "logFile"],
    "QP": ["multiscale", "problemType", "verify", "printLevel", "saveInput", "debug",
           "feasTol", "optTol", "logFile", "solver"],
}


def get_cobra_solver_params_options_for_type(problem_type: Union[str, Sequence[str]]) -> List[str]:
    """Parameter names supported for a problem type (or the union over several).

    MATLAB: paramNames = getCobraSolverParamsOptionsForType(problemType)
    """
    if not isinstance(problem_type, str):              # MATLAB: iscell(problemType)
        names = []
        for t in problem_type:
            names += get_cobra_solver_params_options_for_type(t)
        return sorted(set(names))                       # MATLAB: unique() sorts
    if problem_type not in _OPTIONS_FOR_TYPE:
        raise NotImplementedError(
            f"Parameter list for {problem_type!r} not yet translated (only LP and QP read from source).")
    return list(_OPTIONS_FOR_TYPE[problem_type])


def _default_values(problem_type: str) -> dict:
    """MATLAB valDef: defaults based on Gurobi LP solver defaults."""
    val_def = {
        "feasTol": 1e-6,        # (primal) feasibility tolerance
        "optTol": 1e-6,         # (dual) optimality tolerance
        "objTol": 1e-6,         # only for comparing two objective values
        "minNorm": None,        # MATLAB []
        "printLevel": 0,
        "verify": 0,
        "primalOnly": 0,
        "timeLimit": 1e36,
        "iterationLimit": 1000,
        "logFile": None,        # MATLAB []: empty so no log file is created by default
        "saveInput": None,      # MATLAB []
        "problemType": problem_type,
        "debug": 0,
        "lifting": 0,
        "multiscale": 0,        # true if problem is multiscale
        "intTol": 1e-12,
        "relMipGapTol": 1e-12,
        "absMipGapTol": 1e-12,
        "NUMERICALEMPHASIS": 1,
    }
    if problem_type in ("EP", "CLP"):
        val_def.update(feasTol=1e-6, optTol=1e-6, solver="mosek")
    return val_def


def get_cobra_solver_params(problem_type: str,
                            param_names: Optional[Union[str, Sequence[str]]] = None,
                            param_structure: Optional[Union[dict, str]] = None) -> tuple:
    """Look up solver parameters: param_structure > global > default.

    MATLAB: varargout = getCobraSolverParams(problemType, paramNames, paramStructure)

    Parameters
    ----------
    problem_type : 'LP', 'QP', 'MILP', 'MIQP', 'EP', 'CLP' or 'NLP'
    param_names : name or list of names; None or empty = all names for the type
    param_structure : dict of values to use, or 'default'. If omitted (None),
        only the built-in defaults are used (MATLAB: globals ignored too).

    Returns
    -------
    tuple
        One value per name, followed by the output parameter dict when
        param_structure was given (MATLAB paramStructureOut).
    """
    # --- names (MATLAB: if ~exist('paramNames','var') || isempty(paramNames)) ---
    default_param_names = None
    if not param_names:
        default_param_names = get_cobra_solver_params_options_for_type(problem_type)
        param_names = []
    elif isinstance(param_names, str):
        param_names = [param_names]
    else:
        param_names = list(param_names)

    # --- param structure (MATLAB: if exist('paramStructure','var')) ---
    if param_structure is not None:
        override_defaults = True
        structure = param_structure if isinstance(param_structure, dict) else {}
        param_structure_out = dict(structure)   # MATLAB: paramStructureOut = paramStructure
    else:
        override_defaults = False
        structure = {}
        param_structure_out = {}

    # --- defaults and globals for the problem type ---
    if problem_type not in PROBLEM_TYPES:
        raise ValueError(f"Unrecognized solver type {problem_type!r}")
    val_def = _default_values(problem_type)
    parameters_global = GLOBAL_SOLVER_PARAMS[problem_type]

    # --- no names given: output structure filled from defaults ---
    if default_param_names is not None:
        for name in default_param_names:
            if name in val_def:
                # Quirk kept: overwrites values supplied in param_structure.
                param_structure_out[name] = val_def[name]

    # --- each requested name: default, then global, then supplied value ---
    values = []
    for name in param_names:
        value = None
        if name in val_def:
            value = val_def[name]
            param_structure_out[name] = value
        if override_defaults:
            if name in parameters_global:
                value = parameters_global[name]
                param_structure_out[name] = value
            if name in structure:
                value = structure[name]
                param_structure_out[name] = value
        values.append(value)

    if override_defaults:
        values.append(param_structure_out)
    return tuple(values)
