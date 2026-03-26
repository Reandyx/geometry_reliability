% ========================================
%   MAIN TEST 6 — SORM + GEOMETRY VALIDATION (M6)
% ========================================

clc;
clear;
close all;
addpath(genpath(pwd));
rng(1);

fprintf('\n========================================\n');
fprintf('   MAIN TEST 6 — SORM VALIDATION (M6)\n');
fprintf('========================================\n');

% =========================
% Run study
% =========================
model_type = 'global';

results = run_curvature_study(model_type);

n = length(results);

% =========================
% Extract data
% =========================
a_vals = zeros(n,1);
kappa_global = zeros(n,1);
N_vals = zeros(n,1);
beta = zeros(n,1);

Pf_mc = zeros(n,1);
Pf_form = zeros(n,1);
Pf_sorm = zeros(n,1);

err_form = zeros(n,1);
err_sorm = zeros(n,1);
improvement = zeros(n,1);

for i = 1:n
    a_vals(i) = results(i).a;

    kappa_global(i) = results(i).kappa_global;
    N_vals(i) = results(i).N;

    beta(i) = results(i).beta_FORM;

    Pf_mc(i) = results(i).Pf_MC;
    Pf_form(i) = results(i).Pf_FORM;
    Pf_sorm(i) = results(i).Pf_SORM;

    err_form(i) = results(i).error_FORM;
    err_sorm(i) = results(i).error_SORM;
    improvement(i) = results(i).improvement;
end

% =========================
% GEOMETRY CLASSIFICATION (PER CASE)
% =========================
topology_flag = 0; % M6 = connected geometry

geom_class = strings(n,1);

fprintf('\n--- GEOMETRY CLASSIFICATION DEBUG ---\n');

for i = 1:n
    geom_class(i) = classify_geometry( ...
        N_vals(i), ...
        kappa_global(i), ...
        topology_flag);

    fprintf('Case %d | N = %.3f | kappa = %.3f | class = %s\n', ...
        i, N_vals(i), kappa_global(i), geom_class(i));
end

% =========================
% Sort by normalized curvature
% =========================
[N_vals, idx] = sort(N_vals);

a_vals = a_vals(idx);
kappa_global = kappa_global(idx);
beta = beta(idx);

Pf_mc = Pf_mc(idx);
Pf_form = Pf_form(idx);
Pf_sorm = Pf_sorm(idx);

err_form = err_form(idx);
err_sorm = err_sorm(idx);
improvement = improvement(idx);

geom_class = geom_class(idx); 

% =========================
% PRINT SUMMARY TABLE
% =========================
fprintf('\n========================================\n');
fprintf('   SUMMARY TABLE\n');
fprintf('========================================\n');

fprintf(' a     N         Pf_MC      Pf_FORM    Pf_SORM    err_FORM   err_SORM   class\n');

for i = 1:n
    fprintf('%4.2f  %10.3e  %10.3e  %10.3e  %10.3e  %8.3f   %8.3f   %s\n', ...
        a_vals(i), N_vals(i), ...
        Pf_mc(i), Pf_form(i), Pf_sorm(i), ...
        err_form(i), err_sorm(i), ...
        geom_class(i));
end

% =========================
% VISUALIZATION — LIMIT STATE (M6)
% =========================
fprintf('\n=== VISUALIZATION (M6) ===\n');

% Plot a few representative cases (low, mid, high curvature)
idx_vis = unique([1, round(n/2), n]);

for k = 1:length(idx_vis)
    i = idx_vis(k);

    fprintf('Plotting case a = %.2f (N = %.3f)\n', ...
        results(i).a, results(i).N);

    plot_limit_state_2d(results(i).problem, results(i).U_star);
end

% =========================
% PLOT 1 — ERROR vs CURVATURE
% =========================
figure;
plot(N_vals, err_form, 'o-', 'LineWidth', 2); hold on;
plot(N_vals, err_sorm, 's-', 'LineWidth', 2);

xlabel('Normalized Curvature N = ||H|| / ||\nabla g||');
ylabel('Relative Error');
title('FORM vs SORM Error vs Normalized Curvature');
legend('FORM','SORM','Location','best');
grid on;

% =========================
% PLOT 2 — IMPROVEMENT
% =========================
figure;
plot(N_vals, improvement, 'd-', 'LineWidth', 2);

xlabel('Normalized Curvature N');
ylabel('Error Reduction (FORM - SORM)');
title('SORM Improvement over FORM');
grid on;

% =========================
% PLOT 3 — Pf comparison
% =========================
figure;
semilogy(N_vals, Pf_mc, 'k-o', 'LineWidth', 2); hold on;
semilogy(N_vals, Pf_form, 'r-s', 'LineWidth', 2);
semilogy(N_vals, Pf_sorm, 'b-d', 'LineWidth', 2);

xlabel('Normalized Curvature N');
ylabel('Failure Probability (log scale)');
title('Pf Comparison: MC vs FORM vs SORM');
legend('MC','FORM','SORM','Location','best');
grid on;

% =========================
% PLOT 4 — ERROR vs BETA
% =========================
figure;
plot(beta, err_form, 'o-', 'LineWidth', 2); hold on;
plot(beta, err_sorm, 's-', 'LineWidth', 2);

xlabel('Reliability Index \beta');
ylabel('Relative Error');
title('Error vs Reliability Index');
legend('FORM','SORM','Location','best');
grid on;

% =========================
% QUICK INTERPRETATION
% =========================
fprintf('\n========================================\n');
fprintf('   QUICK INTERPRETATION\n');
fprintf('========================================\n');

for i = 1:n
    
    if err_form(i) < 0.05
        regime = 'LOW curvature (FORM OK)';
    elseif improvement(i) > 0.1
        regime = 'MODERATE curvature (SORM effective)';
    elseif improvement(i) > 0
        regime = 'HIGH curvature (SORM limited)';
    else
        regime = 'FAILURE regime (both unreliable)';
    end
    
    fprintf('a=%.2f | N=%.3e | class=%s | %s\n', ...
        a_vals(i), N_vals(i), geom_class(i), regime);
end

fprintf('\n=== M6 COMPLETE ===\n');