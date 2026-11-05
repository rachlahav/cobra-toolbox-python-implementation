# COBRA Toolbox → Python

Python translation of the COBRA Toolbox whole-body modelling functions
([`src/analysis/wholeBody/PSCMToolbox`](https://github.com/opencobra/cobratoolbox/tree/master/src/analysis/wholeBody/PSCMToolbox))
and the FBA functions they depend on, built on [cobrapy](https://cobrapy.readthedocs.io).
Zomorrodi lab, MGH.

- **Status of every function:** [STATUS.md](STATUS.md)
- **Where cobrapy already does the job:** [docs/cobrapy_equivalents.md](docs/cobrapy_equivalents.md)
- **Latest test findings for `optimize_cb_model`:** [docs/optimize_cb_model_findings.md](docs/optimize_cb_model_findings.md)
- **Basic whole-body model test (Harvey 1.03d):** [docs/wholebody_basic_test_findings.md](docs/wholebody_basic_test_findings.md)

## Layout

```
pscm/                     the translation (package)
  core/                   optimize_cb_model, solver_params, cobra_solver, build_opt_problem
  analysis/  constraints/  data/
tests/
  test_optimize_cb_model.py   checks vs MATLAB reference values (iAF1260, E. coli core)
  adapters.py                 connects implementations to the tests
  reference/                  SciPy/HiGHS port of optimizeCbModel, used as a control
  models/                     ecoli_core_model.mat, iAF1260.mat (from opencobra/COBRA.models)
benchmarks/
  benchmark.py  compare.py    Python timings; side-by-side with MATLAB
  matlab/run_matlab_reference.m
  wholebody/                  basic whole-body test on Harvey 1.03d (MATLAB + Python)
tools/export_workbench.py  sync translations from the PSCM workbench artifact
```

## Run the tests

```bash
pip install -e ".[test]"
pytest -q                                  # translation + reference control
OPTCB_IMPLS=translation pytest -q          # translation only
OPTCB_SOLVER=gurobi pytest -q              # choose the cobrapy solver
```

Every push runs the tests on GitHub Actions (cobrapy with GLPK). The reference
step must always pass: if it fails, the tests or the environment are wrong,
not the translation.

### Where the expected values come from

COBRA's own MATLAB test, `test/verifiedTests/analysis/testOptimizeCbModel/testOptimizeCbModel.m`,
asserts on iAF1260 (tol 1e-6): max growth `0.736700938697735` and weighted L1
total flux `600.3485501480500` (`g1 = 1:n`). The E. coli core model has no
MATLAB values yet; those tests skip until you add numbers from the MATLAB run below.

## Compare with MATLAB

```matlab
% in MATLAB, after initCobraToolbox(false)
cd benchmarks/matlab
run_matlab_reference        % set solverName at the top to match Python
```

```bash
OPTCB_IMPLS=translation python benchmarks/benchmark.py 10 benchmarks/python_results.csv
python benchmarks/compare.py benchmarks/matlab/matlab_results.csv benchmarks/python_results.csv
```

Use the same LP solver on both sides. Plain FBA has alternate optima, so
compare the objective there; the L1 (`minNorm='one'`) flux vectors are the
comparable ones.

## Workflow with the workbench

Translations are written in the PSCM workbench artifact. To bring changes here,
ask Claude to export the workbench's `functions` collection and run
`python tools/export_workbench.py <export dir>`, which rewrites the modules,
`STATUS.md` and `docs/cobrapy_equivalents.md`. Exported files carry a header
naming their workbench record.

## Basic whole-body test (Harvey 1.03d)

Harvey is not stored in this repo (it ships with the COBRA Toolbox at
`test/models/mat/Harvey_1_03d.mat`). Findings so far:
[docs/wholebody_basic_test_findings.md](docs/wholebody_basic_test_findings.md).

```bash
# Python only (numpy + scipy): feasibility, max and min of Whole_body_objective_rxn
python benchmarks/wholebody/wbm_basic_python.py path/to/Harvey_1_03d.mat
python benchmarks/wholebody/diagnose_infeasibility.py path/to/Harvey_1_03d.mat
```

```matlab
% MATLAB vs Python: after initCobraToolbox(false); set solverName at the top
run benchmarks/wholebody/run_matlab_fba.m   % writes harvey_lean.mat, matlab_result.mat
```

```bash
python benchmarks/wholebody/run_python_fba.py glpk   # cobrapy on the same LP, PASS/FAIL vs MATLAB
```
