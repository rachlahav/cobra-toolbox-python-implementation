% run_matlab_fba.m
% Most basic whole-body model (WBM) test, MATLAB side.
%
% Maximises Whole_body_objective_rxn on Harvey 1.03d with optimizeCbModel and
% writes, next to this script:
%   matlab_result.mat  f, stat, origStat, time_s, v  (read by run_python_fba.py)
%   harvey_lean.mat    S, lb, ub, c, b, csense, C, d, dsense, rxns, mets, ctrs
%                      (-v7, so scipy.io.loadmat can read it)
%
% Before running:  initCobraToolbox(false)
% Harvey 1.03d ships with the COBRA Toolbox at test/models/mat/Harvey_1_03d.mat.
%
% Solver: use CPLEX or Gurobi. GLPK on an 8 GB laptop ran for hours without
% finishing (2026-11-02); pdco (interior point) is too loose for a parity test.

solverName = 'ibm_cplex';       % or 'gurobi', 'glpk', ...
modelFile  = fullfile(fileparts(which('initCobraToolbox')), ...
                      'test', 'models', 'mat', 'Harvey_1_03d.mat');
outDir     = fileparts(mfilename('fullpath'));

S = load(modelFile, 'male');
model = S.male; clear S
wbo = 'Whole_body_objective_rxn';
model = changeRxnBounds(model, wbo, 0, 'l');
model = changeRxnBounds(model, wbo, 1000, 'u');
model = changeObjective(model, wbo);
model.osenseStr = 'max';

changeCobraSolver(solverName, 'LP', 0);
t0 = tic; sol = optimizeCbModel(model); t = toc(t0);
fprintf('MATLAB (%s): stat=%d origStat=%s f=%.10g time=%.1fs\n', ...
        solverName, sol.stat, string(sol.origStat), sol.f, t);

lean = struct('S', model.S, 'lb', model.lb, 'ub', model.ub, 'c', model.c, ...
    'b', model.b, 'csense', model.csense, 'C', model.C, 'd', model.d, ...
    'dsense', model.dsense, 'rxns', {model.rxns}, 'mets', {model.mets}, ...
    'ctrs', {model.ctrs});
save(fullfile(outDir, 'harvey_lean.mat'), '-struct', 'lean', '-v7');

res = struct('f', sol.f, 'stat', sol.stat, 'origStat', string(sol.origStat), ...
             'solver', solverName, 'time_s', t, 'v', sol.v);
save(fullfile(outDir, 'matlab_result.mat'), '-struct', 'res', '-v7');
