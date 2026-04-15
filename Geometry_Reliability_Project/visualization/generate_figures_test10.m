clear; clc;

addpath(genpath(fileparts(pwd)));

assert(exist(fullfile('results', 'results_main_test10.mat'), 'file') == 2, ...
    'Run main_test10 first.');

load('results/results_main_test10.mat');

%% ========================================
% SORT DATA (CRITICAL)
% =========================================

% --- CURVATURE ---
[gamma_vals, idx] = sort([results.curvature.gamma]);
err_FORM = [results.curvature.err_form];
Pf_FORM_c = [results.curvature.Pf_form];
Pf_MCS_c  = [results.curvature.Pf_mcs];
Pf_SORM   = [results.curvature.Pf_sorm];
kappa_vals = [results.curvature.kappa];

err_FORM = err_FORM(idx);
Pf_FORM_c = Pf_FORM_c(idx);
Pf_MCS_c  = Pf_MCS_c(idx);
Pf_SORM   = Pf_SORM(idx);
kappa_vals = kappa_vals(idx);

% --- TOPOLOGY ---
[c_vals, idx] = sort([results.topology.c]);
Pf_FORM = [results.topology.Pf_form];
Pf_MCS  = [results.topology.Pf_mcs];
err_topo = [results.topology.error];

Pf_FORM = Pf_FORM(idx);
Pf_MCS  = Pf_MCS(idx);
err_topo = err_topo(idx);

% --- RARE EVENTS ---
[beta_rare, idx] = sort([results.rare.beta]);
N_MCS = [results.rare.N_mcs];
N_IS  = [results.rare.N_is];
CoV_IS = [results.rare.CoV_is];

N_MCS = N_MCS(idx);
N_IS  = N_IS(idx);
CoV_IS = CoV_IS(idx);

%% ========================================
% FIGURE 1 — CURVATURE EFFECT
% =========================================
figure;
scatter(gamma_vals, err_FORM, 70, 'filled'); hold on;

p = polyfit(gamma_vals, err_FORM, 2);
x_fit = linspace(min(gamma_vals), max(gamma_vals), 200);
y_fit = polyval(p, x_fit);
plot(x_fit, y_fit, '--', 'LineWidth', 1.5);

grid on;
xlabel('\gamma = \beta \kappa');
ylabel('Relative FORM Error');
title('FORM Error vs Curvature');

%% ========================================
% FIGURE 2 — GEOMETRY CLASSIFICATION
% =========================================
gamma_vals = [results.curvature.gamma];
err_FORM   = [results.curvature.err_form];   

figure;
scatter(gamma_vals, err_FORM, 70, 'filled'); hold on;

xlim([0 max(gamma_vals)*1.075]);

xline(0.2,'--','LineWidth',1.2);
xline(0.6,'--','LineWidth',1.2);

grid on;
xlabel('\gamma = \beta \kappa');
ylabel('Relative FORM Error');
title('Geometry-Based Classification');

text(0.05, 0.03, 'Class I (<5%)');
text(0.25, 0.15, 'Class II (5–30%)');
text(1.0, 0.2, 'Class III (>30%)');

%% ========================================
% FIGURE 3 — TOPOLOGY FAILURE (Pf)
% =========================================
figure;
semilogy(c_vals, Pf_FORM, '-o', 'LineWidth', 2); hold on;
semilogy(c_vals, Pf_MCS, '-s', 'LineWidth', 2);

grid on;
xlabel('Topology parameter c');
ylabel('Failure Probability P_f');
title('FORM vs MCS — Topological Failure');

legend('FORM', 'MCS', 'Location', 'southwest');

%% ========================================
% FIGURE 4 — TOPOLOGY ERROR
% =========================================
figure;
plot(c_vals, err_topo, '-o', 'LineWidth', 2);

grid on;
xlabel('Topology parameter c');
ylabel('Relative Error');
title('FORM Error due to Topology');

%% ========================================
% FIGURE 5 — RARE EVENT SCALING
% =========================================
figure;
semilogy(beta_rare, N_MCS, '-o', 'LineWidth', 2); hold on;
semilogy(beta_rare, N_IS, '-s', 'LineWidth', 2);

grid on;
xlabel('\beta');
ylabel('Required Samples');
title('Sample Complexity: MCS vs IS');

legend('MCS', 'Importance Sampling', 'Location', 'northwest');

%% ========================================
% FIGURE 6 — IS STABILITY
% =========================================
figure;
plot(beta_rare, CoV_IS, '-o', 'LineWidth', 2);

grid on;
xlabel('\beta');
ylabel('Coefficient of Variation');
title('Importance Sampling Stability');

%% ========================================
% FIGURE 7 — CURVATURE VALIDATION
% =========================================
figure;
plot(gamma_vals, kappa_vals, '-o', 'LineWidth', 2);

grid on;
xlabel('\gamma');
ylabel('Principal Curvature \kappa');
title('Curvature at Design Point');

%% ========================================
% FIGURE 8 — FORM vs SORM vs MCS
% =========================================
figure;
semilogy(gamma_vals, Pf_FORM_c, '-o', 'LineWidth', 2); hold on;
semilogy(gamma_vals, Pf_SORM, '-s', 'LineWidth', 2);
semilogy(gamma_vals, Pf_MCS_c, '-d', 'LineWidth', 2);

grid on;
xlabel('\gamma');
ylabel('Failure Probability P_f');
title('FORM vs SORM vs MCS');

legend('FORM','SORM','MCS','Location','southwest');

%% ========================================
% FIGURE 9 — TOPOLOGY FAILURE RATIO (CRITICAL)
% =========================================
ratio = Pf_FORM ./ Pf_MCS;

figure;
semilogy(c_vals, ratio, '-o', 'LineWidth', 2);

grid on;
xlabel('Topology parameter c');
ylabel('P_f^{FORM} / P_f^{MCS}');
title('Topology-Induced Failure Bias');