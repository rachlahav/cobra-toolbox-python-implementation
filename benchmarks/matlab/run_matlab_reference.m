% run_matlab_reference.m
% Runs optimizeCbModel on the same cases as benchmark.py and writes
% matlab_results.csv (same columns) plus matlab_results_fluxes.mat.
%
% Before running:  initCobraToolbox(false); then cd into this folder.
% Set the solver below to the one you will also use from Python if
% possible (e.g. 'gurobi' on both sides) so timings are comparable.

solverName = 'gurobi';      % or 'glpk', 'ibm_cplex', 'mosek', ...
repeats    = 10;
modelDir   = fullfile('..', '..', 'tests', 'models');

changeCobraSolver(solverName, 'LP', 0);

models = {'ecoli_core', 'ecoli_core_model.mat', 'model';
          'iAF1260',    'iAF1260.mat',          'iAF1260'};
cases  = {'fba_max',     'max', 0,     false;
          'l1_one',      'max', 'one', false;
          'l1_weighted', 'max', 'one', true};

fid = fopen('matlab_results.csv', 'w');
fprintf(fid, 'model,case,stat,f,sum_abs_v,mean_s,sd_s,median_s,repeats,load_s,impl\n');
fluxes = struct();

for i = 1:size(models, 1)
    t0 = tic;
    S = load(fullfile(modelDir, models{i, 2}));
    base = S.(models{i, 3});
    loadS = toc(t0);
    for k = 1:size(cases, 1)
        m = base;
        if cases{k, 4}
            m.g1 = (1:size(m.S, 2))';
        end
        optimizeCbModel(m, cases{k, 2}, cases{k, 3});   % warm-up, not timed
        t = zeros(repeats, 1);
        for r = 1:repeats
            t1 = tic;
            sol = optimizeCbModel(m, cases{k, 2}, cases{k, 3});
            t(r) = toc(t1);
        end
        fprintf(fid, '%s,%s,%d,%.15g,%.15g,%.6f,%.6f,%.6f,%d,%.4f,%s\n', ...
            models{i, 1}, cases{k, 1}, sol.stat, sol.f, sum(abs(sol.v)), ...
            mean(t), std(t), median(t), repeats, loadS, ['matlab-' solverName]);
        fluxes.([models{i, 1} '__' cases{k, 1}]) = sol.v;
        fprintf('%-10s %-12s f=%.10f sum|v|=%.6f median %8.2f ms\n', ...
            models{i, 1}, cases{k, 1}, sol.f, sum(abs(sol.v)), 1e3 * median(t));
    end
end
fclose(fid);
save('matlab_results_fluxes.mat', '-struct', 'fluxes');
fprintf('\nwrote matlab_results.csv and matlab_results_fluxes.mat\n');
