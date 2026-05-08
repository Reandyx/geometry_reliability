clc;
clear;
close all;
addpath(genpath(pwd));
rng(1);

fprintf('\n========================================\n');
fprintf('   MAIN TEST 5 — GEOMETRY STUDY (M5)\n');
fprintf('========================================\n');

results_global = run_curvature_study('global');
results_local  = run_curvature_study('local');

T_global = struct2table(results_global);
T_local  = struct2table(results_local);

fprintf('\n============================\n');
fprintf('   GLOBAL RESULTS\n');
fprintf('============================\n');
disp(T_global);

fprintf('\n============================\n');
fprintf('   LOCAL RESULTS\n');
fprintf('============================\n');
disp(T_local);

%% =========================
% M4 PLOTS (REUSED FROM M4)
% These plots are intentionally kept here for comparison/reporting.
% =========================

figure; hold on; grid on;
plot([results_global.a], [results_global.error_FORM], '-o', 'LineWidth', 2);
plot([results_local.a],  [results_local.error_FORM],  '-s', 'LineWidth', 2);
xlabel('a');
ylabel('Relative Error');
legend('Global', 'Local', 'Location', 'best');
title('FORM Error vs Curvature Parameter');

figure; hold on; grid on;
plot([results_global.a], [results_global.bias], '-o', 'LineWidth', 2);
plot([results_local.a],  [results_local.bias],  '-s', 'LineWidth', 2);
xlabel('a');
ylabel('Bias');
legend('Global', 'Local', 'Location', 'best');
title('FORM Bias vs Curvature Parameter');

figure; hold on; grid on;
plot([results_global.a], [results_global.beta_FORM], '-o', 'LineWidth', 2);
plot([results_local.a],  [results_local.beta_FORM],  '-s', 'LineWidth', 2);
xlabel('a');
ylabel('\beta');
legend('Global', 'Local', 'Location', 'best');
title('Reliability Index vs Curvature Parameter');

figure; hold on; grid on;
errorbar([results_global.a], [results_global.Pf_MC], [results_global.sigma_MC], '-o', 'LineWidth', 2);
plot([results_global.a], [results_global.Pf_FORM], '-o', 'LineWidth', 2);

errorbar([results_local.a], [results_local.Pf_MC], [results_local.sigma_MC], '-s', 'LineWidth', 2);
plot([results_local.a], [results_local.Pf_FORM], '-s', 'LineWidth', 2);

xlabel('a');
ylabel('Failure Probability');
legend('MC Global', 'FORM Global', 'MC Local', 'FORM Local', 'Location', 'best');
title('Pf Comparison');

figure; hold on; grid on;
U_global = reshape([results_global.U_star], 2, []);
U_local  = reshape([results_local.U_star],  2, []);
plot(U_global(1,:), U_global(2,:), '-o', 'LineWidth', 2);
plot(U_local(1,:),  U_local(2,:),  '-s', 'LineWidth', 2);
xlabel('U_1^*');
ylabel('U_2^*');
legend('Global', 'Local', 'Location', 'best');
title('Design Point Trajectory');

%% =========================
% M5 NEW ANALYSIS
% =========================

% Required fields for M5 main analysis
assert(isfield(results_global, 'kappa_global'), ...
    'run_curvature_study(''global'') must return field: kappa_global');
assert(isfield(results_local, 'kappa_global'), ...
    'run_curvature_study(''local'') must return field: kappa_global');
assert(isfield(results_global, 'error_FORM'), ...
    'run_curvature_study(''global'') must return field: error_FORM');
assert(isfield(results_local, 'error_FORM'), ...
    'run_curvature_study(''local'') must return field: error_FORM');

% Extract main M5 data
kappa_g = [results_global.kappa_global];
error_g = [results_global.error_FORM];

kappa_l = [results_local.kappa_global];
error_l = [results_local.error_FORM];

% Remove invalid values
valid_g = isfinite(kappa_g) & isfinite(error_g) & (kappa_g > 0) & (error_g > 0);
valid_l = isfinite(kappa_l) & isfinite(error_l) & (kappa_l > 0) & (error_l > 0);

kappa_g_plot = kappa_g(valid_g);
error_g_plot = error_g(valid_g);

kappa_l_plot = kappa_l(valid_l);
error_l_plot = error_l(valid_l);

%% =========================
% MAIN FIGURE — Error vs kappa
% =========================
figure; hold on; grid on;

scatter(kappa_g_plot, error_g_plot, 80, 'o', 'filled');
scatter(kappa_l_plot, error_l_plot, 80, 's', 'filled');

set(gca, 'XScale', 'log', 'YScale', 'log');

xlabel('\kappa (curvature)');
ylabel('FORM relative error');
legend('Global', 'Local', 'Location', 'best');
title('FORM Error vs Curvature (M5 MAIN RESULT)');

%% =========================
% SECONDARY — Error vs N
% Plot only if field N exists in both result sets
% =========================
hasN_global = isfield(results_global, 'N');
hasN_local  = isfield(results_local, 'N');

if hasN_global && hasN_local
    N_g = [results_global.N];
    N_l = [results_local.N];

    validNg = valid_g & isfinite(N_g) & (N_g > 0);
    validNl = valid_l & isfinite(N_l) & (N_l > 0);

    N_g_plot = N_g(validNg);
    err_g_N  = error_g(validNg);

    N_l_plot = N_l(validNl);
    err_l_N  = error_l(validNl);

    figure; hold on; grid on;

    scatter(N_g_plot, err_g_N, 80, 'o', 'filled');
    scatter(N_l_plot, err_l_N, 80, 's', 'filled');

    set(gca, 'XScale', 'log', 'YScale', 'log');

    xlabel('Nonlinearity Index N');
    ylabel('FORM relative error');
    legend('Global', 'Local', 'Location', 'best');
    title('FORM Error vs Nonlinearity Index');
else
    warning(['Skipping "Error vs N" plot: field N is missing in results. ', ...
             'Add field N in run_curvature_study() if this plot is required.']);
end

%% =========================
% OPTIONAL — Curvature diagnostics
% Plot principal curvature only if available
% =========================
hasKP_global = isfield(results_global, 'kappa_principal');
hasKP_local  = isfield(results_local, 'kappa_principal');

figure; hold on; grid on;

% Global total/global curvature
plot([results_global.a], [results_global.kappa_global], '-o', 'LineWidth', 2);

% Local total/global curvature
plot([results_local.a], [results_local.kappa_global], '-s', 'LineWidth', 2);

legendEntries = {'Global |κ|', 'Local |κ|'};

% Principal curvature only if field exists
if hasKP_global
    plot([results_global.a], [results_global.kappa_principal], '--o', 'LineWidth', 2);
    legendEntries{end+1} = 'Global \kappa_{principal}';
end

if hasKP_local
    plot([results_local.a], [results_local.kappa_principal], '--s', 'LineWidth', 2);
    legendEntries{end+1} = 'Local \kappa_{principal}';
end

xlabel('a');
ylabel('\kappa');
legend(legendEntries, 'Location', 'best');
title('Curvature Evolution vs Parameter a');

fprintf('\n========================================\n');
fprintf('   MAIN TEST 5 COMPLETE\n');
fprintf('========================================\n');