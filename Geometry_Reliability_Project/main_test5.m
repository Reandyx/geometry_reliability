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
% M4 PLOTS (UNCHANGED)
% =========================

figure; hold on; grid on;
plot([results_global.a], [results_global.error_FORM], '-o', 'LineWidth', 2);
plot([results_local.a],  [results_local.error_FORM],  '-s', 'LineWidth', 2);
xlabel('a');
ylabel('Relative Error');
legend('Global', 'Local');
title('FORM Error vs Curvature');

figure; hold on; grid on;
plot([results_global.a], [results_global.bias], '-o', 'LineWidth', 2);
plot([results_local.a],  [results_local.bias],  '-s', 'LineWidth', 2);
xlabel('a');
ylabel('Bias');
legend('Global', 'Local');
title('FORM Bias vs Curvature');

figure; hold on; grid on;
plot([results_global.a], [results_global.beta_FORM], '-o', 'LineWidth', 2);
plot([results_local.a],  [results_local.beta_FORM],  '-s', 'LineWidth', 2);
xlabel('a');
ylabel('\beta');
legend('Global', 'Local');
title('Reliability Index vs Curvature');

figure; hold on; grid on;
errorbar([results_global.a], [results_global.Pf_MC], [results_global.sigma_MC], '-o', 'LineWidth', 2);
plot([results_global.a], [results_global.Pf_FORM], '-o', 'LineWidth', 2);

errorbar([results_local.a], [results_local.Pf_MC], [results_local.sigma_MC], '-s', 'LineWidth', 2);
plot([results_local.a], [results_local.Pf_FORM], '-s', 'LineWidth', 2);

xlabel('a');
ylabel('Failure Probability');
legend('MC Global','FORM Global','MC Local','FORM Local');
title('Pf Comparison');

figure; hold on; grid on;
U_global = reshape([results_global.U_star], 2, []);
U_local  = reshape([results_local.U_star],  2, []);
plot(U_global(1,:), U_global(2,:), '-o', 'LineWidth', 2);
plot(U_local(1,:),  U_local(2,:),  '-s', 'LineWidth', 2);
xlabel('U_1^*');
ylabel('U_2^*');
legend('Global','Local');
title('Design Point Trajectory');

%% =========================
%  M5 NEW PLOTS (FINAL)
% =========================

% Extract data
kappa_g = [results_global.kappa_global];
error_g = [results_global.error_FORM];
N_g     = [results_global.N];

kappa_l = [results_local.kappa_global];
error_l = [results_local.error_FORM];
N_l     = [results_local.N];

% Remove invalid values
valid_g = isfinite(kappa_g) & isfinite(error_g) & error_g > 0;
valid_l = isfinite(kappa_l) & isfinite(error_l) & error_l > 0;

kappa_g = kappa_g(valid_g);
error_g = error_g(valid_g);
N_g = N_g(valid_g);

kappa_l = kappa_l(valid_l);
error_l = error_l(valid_l);
N_l = N_l(valid_l);

%% =========================
% MAIN FIGURE — Error vs kappa
% =========================
figure; hold on; grid on;

scatter(kappa_g, error_g, 80, 'o', 'filled');
scatter(kappa_l, error_l, 80, 's', 'filled');

set(gca, 'XScale', 'log', 'YScale', 'log');

xlabel('\kappa (curvature)');
ylabel('FORM relative error');
legend('Global','Local','Location','best');

title('FORM Error vs Curvature (M5 MAIN RESULT)');

%% =========================
% SECONDARY — Error vs N
% =========================
figure; hold on; grid on;

scatter(N_g, error_g, 80, 'o', 'filled');
scatter(N_l, error_l, 80, 's', 'filled');

set(gca, 'XScale', 'log', 'YScale', 'log');

xlabel('Nonlinearity Index N');
ylabel('FORM relative error');
legend('Global','Local','Location','best');

title('FORM Error vs Nonlinearity Index');

%% =========================
% OPTIONAL — kappa vs a (diagnostic)
% =========================
figure; hold on; grid on;

% Global
plot([results_global.a], [results_global.kappa_global], '-o', 'LineWidth', 2);
plot([results_global.a], [results_global.kappa_principal], '--o', 'LineWidth', 2);

% Local
plot([results_local.a], [results_local.kappa_global], '-s', 'LineWidth', 2);
plot([results_local.a], [results_local.kappa_principal], '--s', 'LineWidth', 2);

xlabel('a');
ylabel('\kappa');
legend('Global |κ|','Global κ','Local |κ|','Local κ');

title('Curvature Evolution vs Parameter a');