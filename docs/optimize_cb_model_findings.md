# optimize_cb_model: test findings

Translation under test: `pscm/core/optimize_cb_model.py` (workbench record
`optimizeCbModel`, version 3, updated 2026-10-12).

## 1. Crashes on solvers without an `optimality` tolerance, including GLPK

`_apply_params` reads and sets `cfg.tolerances.feasibility` and
`cfg.tolerances.optimality` unconditionally. optlang raises `AttributeError`
for tolerance names a solver doesn't define:

| optlang backend | feasibility | optimality |
|---|---|---|
| gurobi, cplex, coinor_cbc | yes | yes |
| **glpk** (cobrapy's default) | yes | **no** |
| scipy | no | no |
| highs | yes | no (`ipm_optimality`, `mip_gap` instead) |

So every call fails before solving unless the model uses Gurobi, CPLEX or CBC.
Observed directly with the SciPy backend; the GLPK case follows from optlang's
source and is what the "Translation (pscm, GLPK)" CI step checks.

Suggested fix: set only the tolerances the solver has, e.g.

```python
tol = cfg.tolerances
available = set(dir(tol))          # SolverTolerances.__dir__ lists the solver's names
old = {k: getattr(tol, k) for k in ("feasibility", "optimality") if k in available}
if "feasibility" in available:
    tol.feasibility = feas_tol
if "optimality" in available:
    tol.optimality = opt_tol
```

and restore only the keys in `old`.

## 2. Slacks need `constraint.primal` (not available on the SciPy backend)

`_solve` computes `slack = rhs - con.primal` for every row. optlang's SciPy
backend can't return constraint primals (`TypeError` inside optlang). This is
an optlang limitation, not a translation bug, but it means the translation
can't run on SciPy. If that matters, compute the slack from the fluxes
instead (`b - S @ v`), which also avoids one solver call per row.

## 3. Behaviour the tests pin down (to check once 1 is fixed)

- MATLAB reference values on iAF1260 (from COBRA's `testOptimizeCbModel.m`):
  growth 0.736700938697735; weighted L1 total flux 600.34855014805.
- `f1` for `minNorm='one'` is the **unweighted** `sum|v|`, even with `g1`.
- Infeasible: `stat = 0`, `f = NaN`.
- The cobrapy model must be unchanged after the call (the test adapter compares
  bounds, objective, direction, variable/constraint counts, solver and notes).
- Duals in COBRA's convention: `osense*c - S'y - w = 0`.

## Environment notes

- Large models are slow in cobrapy without `symengine` (iAF1260: ~7 min to load,
  ~1 min per `model.copy()` in a sandbox without it). `pip install symengine`.
- The test adapter reuses one loaded model per test file and undoes bound
  changes with `with model:`, instead of copying it per call.
