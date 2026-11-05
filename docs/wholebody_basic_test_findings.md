# Basic whole-body model test: findings

The simplest whole-body check: FBA on **Harvey 1.03d** (male WBM shipped with the
COBRA Toolbox at `test/models/mat/Harvey_1_03d.mat`), objective
`Whole_body_objective_rxn`. Scripts are in [`benchmarks/wholebody/`](../benchmarks/wholebody).

Status (2026-11-05): **Harvey 1.03d is feasible as shipped.** With HiGHS's
default settings, the Python LP finds a flux with `Whole_body_objective_rxn = 1` that
satisfies mass balance, all 103,600 coupling constraints and the bounds (max
violations 1.6e-9, 1.0e-7, 7.1e-8). The "infeasible" results from 2026-11-03 were
solver artefacts: tightened tolerances (1e-9) and a HiGHS presolve error, not the
model or the LP setup. The MATLAB side has not produced a result yet, so there is
no MATLAB-vs-Python comparison yet.

Test: [`tests/test_wholebody_basic.py`](../tests/test_wholebody_basic.py) (skipped unless
`HARVEY_MAT` points at the model).

## The model as loaded in Python

| | |
|---|---|
| Metabolites × reactions | 56,452 × 81,094 (nnz S = 283,529) |
| Coupling constraints `C v ≤/≥ d` | 103,600: 70,009 `L`, 33,591 `G` (nnz C = 209,602), all `d = 0` |
| `csense` | all `E`, `b = 0` |
| `Whole_body_objective_rxn` as shipped | lb = ub = 1, `c` all zero, `osenseStr = 'min'` |
| Forced fluxes | 23 reactions with lb > 0, 144 with ub < 0 |
| Load time (`scipy.io.loadmat`) | ~4 s |

## Results

### 2026-11-05: `tests/test_wholebody_basic.py` (SciPy 1.18, `method="highs"`, default tolerances)

| Test | Result |
|---|---|
| Dimensions, senses, objective bounds, data consistency (4 tests) | pass |
| Feasibility, model as shipped (WBO fixed = 1) | **pass**: feasible, WBO = 1.0, ~175 s |
| Relaxed problem (WBO in [0, 1000]) is feasible, checked without a solver | **pass**: the as-shipped flux satisfies it |
| Max WBO, bounds [0, 1000] | xfail: HiGHS reports "infeasible" in < 1 s |
| Min WBO, bounds [0, 1000] | xfail: HiGHS reports "infeasible" in < 1 s |

The relaxed problem contains every solution of the as-shipped one, so the "infeasible"
answers for max/min are wrong. They come back in under a second, i.e. from HiGHS
presolve. With presolve off, the max solve had not finished after ~7 min (stopped).

### 2026-11-03: `benchmarks/wholebody/wbm_basic_python.py` (superseded)

All three variants (as shipped, max, min; dual simplex and interior point) were reported
infeasible in < 0.5 s. That script sets feasibility/optimality tolerances to 1e-9;
with the default tolerances the as-shipped case is feasible (above). The script now uses the defaults. Raw output:
[`benchmarks/wholebody/results/2026-11-03_harvey_1_03d_highs.log`](../benchmarks/wholebody/results/2026-11-03_harvey_1_03d_highs.log).

`diagnose_infeasibility.py` found no data problems: no lb > ub, no NaN/Inf, no row of
`S v = 0` or `C v ≤/≥ 0` impossible on its own given the bounds. Its "infeasible with only
the L / only the G coupling rows" results use the relaxed bounds, so they hit the same
presolve error.

## Why it matters for the translation

The model is badly scaled: bounds reach ±1e6 (lb up to 62,169, ub down to −132,714) while
coupling coefficients are small. Solver settings change the answer:
- tolerances tighter than the default (1e-9) turn a feasible model "infeasible";
- HiGHS presolve declares the relaxed model infeasible.

So parity tests on whole-body models must use the same solver and tolerances on both
sides, and should check feasibility of the returned flux directly (as the test does)
rather than trusting the status alone. COBRA's default `feasTol` is 1e-6.

## Next steps

1. Run [`run_matlab_fba.m`](../benchmarks/wholebody/run_matlab_fba.m) with CPLEX or Gurobi
   to get MATLAB's max `Whole_body_objective_rxn`, then compare with
   [`run_python_fba.py`](../benchmarks/wholebody/run_python_fba.py) (cobrapy, same LP).
2. Get a Python max/min value: try HiGHS with presolve off and a longer time limit, or
   Gurobi/CPLEX from Python.

`wbm_basic_python.py` now uses HiGHS's default tolerances (changed 2026-11-05).

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
