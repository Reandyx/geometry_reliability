clear; clc;

addpath(genpath(fileparts(pwd))); % ensures root access

load(fullfile('results', 'results_verification_suite_main_test9.mat'));

assert(exist(fullfile('results', 'results_verification_suite_main_test9.mat'), 'file') == 2, ...
    'Results file not found. Run main_test9 first.');

fprintf('========================================\n');
fprintf(' GENERATING PLOTS\n');
fprintf('========================================\n');

%% ============================================================
% PLOT 1 — FORM ERROR vs GEOMETRY INDEX γ
% ============================================================
figure;

gamma = arrayfun(@(s) s.gamma, results.geometry);
error = arrayfun(@(s) s.rel_error, results.geometry);

plot(gamma, error, '-o', 'LineWidth', 2);
grid on;

xlabel('\gamma = \beta max|\kappa_i|');
ylabel('Relative Error (FORM)');
title('FORM Error vs Geometry Index');

%% ============================================================
% PLOT 2 — METHOD COMPARISON (FORM vs SORM vs MCS)
% ============================================================
figure;

a_vals = arrayfun(@(s) s.a, results.curvature);
Pf_form = arrayfun(@(s) s.Pf_form, results.curvature);
Pf_sorm = arrayfun(@(s) s.Pf_sorm, results.curvature);
Pf_mcs  = arrayfun(@(s) s.Pf_mcs, results.curvature);

plot(a_vals, Pf_form, '-o', 'LineWidth', 2); hold on;
plot(a_vals, Pf_sorm, '-s', 'LineWidth', 2);
plot(a_vals, Pf_mcs,  '-^', 'LineWidth', 2);

grid on;
legend('FORM', 'SORM', 'MCS', 'Location', 'best');

xlabel('Curvature Parameter a');
ylabel('Failure Probability P_f');
title('Method Comparison Across Curvature');

%% ============================================================
% PLOT 3 — RELATIVE ERROR (FORM vs SORM)
% ============================================================
figure;

err_form = abs(Pf_form - Pf_mcs) ./ Pf_mcs;
err_sorm = abs(Pf_sorm - Pf_mcs) ./ Pf_mcs;

plot(a_vals, err_form, '-o', 'LineWidth', 2); hold on;
plot(a_vals, err_sorm, '-s', 'LineWidth', 2);

grid on;
legend('FORM Error', 'SORM Error', 'Location', 'best');

xlabel('Curvature Parameter a');
ylabel('Relative Error');
title('FORM vs SORM Error');

%% ============================================================
% PLOT 4 — TOPOLOGY COMPARISON (BAR)
% ============================================================
figure;

Pf_topo = [results.topology.Pf_form, results.topology.Pf_mcs];

bar(Pf_topo);

set(gca, 'XTickLabel', {'FORM', 'MCS'});

ylabel('Failure Probability');
title('Disconnected Domain: FORM vs MCS');

grid on;

%% ============================================================
% PLOT 5 — RARE EVENT SCALING
% ============================================================
figure;

beta_vals = arrayfun(@(s) s.beta_target, results.rare_event);
Pf_form = arrayfun(@(s) s.Pf_form, results.rare_event);
Pf_mcs  = arrayfun(@(s) s.Pf_mcs, results.rare_event);
Pf_is   = arrayfun(@(s) s.Pf_is, results.rare_event);

semilogy(beta_vals, Pf_form, '-o', 'LineWidth', 2); hold on;
semilogy(beta_vals, Pf_mcs,  '-s', 'LineWidth', 2);
semilogy(beta_vals, Pf_is,   '-^', 'LineWidth', 2);

grid on;
legend('FORM', 'MCS', 'IS', 'Location', 'best');

xlabel('\beta');
ylabel('Failure Probability (log scale)');
title('Rare Event Estimation');

%% ============================================================
% PLOT 6 — IS PERFORMANCE (CoV & ESS)
% ============================================================
figure;

CoV = arrayfun(@(s) s.CoV_is, results.rare_event);
ESS = arrayfun(@(s) s.ESS, results.rare_event);

yyaxis left
plot(beta_vals, CoV, '-o', 'LineWidth', 2);
ylabel('CoV (IS)');

yyaxis right
plot(beta_vals, ESS, '-s', 'LineWidth', 2);
ylabel('Effective Sample Size');

xlabel('\beta');
title('Importance Sampling Performance');

grid on;

fprintf('\nAll plots generated.\n');

%% ============================================================
% PLOT 7 — REGIME MAP (γ vs ERROR)
% ============================================================
figure;

gamma = arrayfun(@(s) s.gamma, results.geometry);
error = arrayfun(@(s) s.rel_error, results.geometry);

scatter(gamma, error, 80, 'filled');
grid on;

xlabel('\gamma = \beta max|\kappa_i|');
ylabel('FORM Relative Error');
title('Geometry Regime Map');

% Threshold lines (you will refine later)
hold on;
xline(0.1, '--', 'Low curvature');
xline(1.0, '--', 'High curvature');

yline(0.05, '--', '5% error');
yline(0.30, '--', '30% error');