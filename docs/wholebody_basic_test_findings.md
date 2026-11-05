# Basic whole-body model test: findings

The simplest whole-body check: FBA on **Harvey 1.03d** (male WBM shipped with the
COBRA Toolbox at `test/models/mat/Harvey_1_03d.mat`), objective
`Whole_body_objective_rxn`. Scripts are in [`benchmarks/wholebody/`](../benchmarks/wholebody).

Status (2026-11-05): **the Python side runs; the solver reports the model
infeasible. The MATLAB side has not produced a result yet**, so there is no
MATLAB-vs-Python comparison yet.

## The model as loaded in Python

| | |
|---|---|
| Metabolites × reactions | 56,452 × 81,094 (nnz S = 283,529) |
| Coupling constraints `C v ≤/≥ d` | 103,600: 70,009 `L`, 33,591 `G` (nnz C = 209,602), all `d = 0` |
| `csense` | all `E`, `b = 0` |
| `Whole_body_objective_rxn` as shipped | lb = ub = 1, `c` all zero, `osenseStr = 'min'` |
| Forced fluxes | 23 reactions with lb > 0, 144 with ub < 0 |
| Load time (`scipy.io.loadmat`) | ~4 s |

## Results (`wbm_basic_python.py`, SciPy 1.18 HiGHS, 2026-11-03)

| Test | HiGHS dual simplex | HiGHS interior point |
|---|---|---|
| 1. Feasibility, model as shipped (WBO fixed = 1) | infeasible (0.4 s) | infeasible (0.3 s) |
| 2. Max WBO, bounds relaxed to [0, 1000] | infeasible (0.4 s) | infeasible (0.4 s) |
| 3. Min WBO, bounds relaxed to [0, 1000] | infeasible (0.4 s) | infeasible (0.4 s) |

Raw output: [`benchmarks/wholebody/results/2026-11-03_harvey_1_03d_highs.log`](../benchmarks/wholebody/results/2026-11-03_harvey_1_03d_highs.log).
Same result with default tolerances and with feasibility/optimality tolerances of 1e-9.

## Diagnostics (`diagnose_infeasibility.py`)

- No lb > ub, no NaN/Inf, no row of `S v = 0` or `C v ≤/≥ 0` that is impossible on its own
  given the bounds (interval arithmetic).
- Relaxing WBO to 0 does not help, so the conflict is independent of the objective.
- Infeasible with only the `L` coupling rows, and with only the `G` rows.
- Without any coupling rows, HiGHS did not detect infeasibility in presolve (the solve was
  stopped after 10 min, unfinished).

So the conflict most likely involves the coupling constraints together with the forced
fluxes.

## Possible causes (not yet distinguished)

1. **The model file.** This Harvey copy may not be feasible as stored (e.g. it expects a diet
   or other setup before solving).
2. **Numerical scaling.** Bounds reach 1e6 while coupling coefficients are small; HiGHS may
   declare infeasible what a MATLAB solver with looser tolerances accepts.
3. **A mistake in the Python LP setup**, e.g. the sign convention of the coupling rows. Row
   counts and senses match the file; signs have not been checked against MATLAB.

## What settles it

Run [`run_matlab_fba.m`](../benchmarks/wholebody/run_matlab_fba.m) with CPLEX or Gurobi. If
MATLAB also returns `stat = 0`, it is the model file; if MATLAB finds a solution, compare with
[`run_python_fba.py`](../benchmarks/wholebody/run_python_fba.py) (cobrapy, same LP from
MATLAB's exported `harvey_lean.mat`).

Solver notes from the attempt on 2026-11-02 (MacBook Air, Apple Silicon, 8 GB RAM, MATLAB R2026b):
- GLPK ran for hours without finishing and used ~21 GB swap.
- pdco (the default COBRA picked) is an interior-point solver too loose for parity checks.
- CPLEX was configured (`ILOG_CPLEX_PATH = ~/Applications/CPLEX_Studio2212/cplex/matlab/x86-64_osx`)
  but not detected: an `x86-64_osx` (Intel) interface can't load in native Apple Silicon
  MATLAB (`maca64`). Check whether an Apple Silicon interface folder exists, or use Gurobi
  (free academic licence, native arm64).

## For comparison: iAF1260 (non-whole-body) in Python, 2026-11-05

`tests/reference/reference_optimize.py` (SciPy/HiGHS port of `optimizeCbModel`) against the
values COBRA's MATLAB `testOptimizeCbModel.m` asserts:

| Case | Python | MATLAB reference | abs. diff |
|---|---|---|---|
| FBA max growth `f` | 0.7367009388433972 | 0.736700938697735 | 1.5e-10 |
| `minNorm='one'`, `g1 = 1:n`: `sum\|v\|` | 600.3485498460217 | 600.3485501480500 | 3.0e-7 |
| `minNorm='one'`, unweighted: `sum\|v\|` | 599.7682158842495 | 599.768216071444 | 1.9e-7 |

All within the MATLAB test's tolerance (1e-6). This is the reference port, not
`pscm/core/optimize_cb_model.py`; see [optimize_cb_model_findings.md](optimize_cb_model_findings.md)
for the translation.
