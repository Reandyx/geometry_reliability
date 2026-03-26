% ========================================
%   MAIN TEST 7 — TOPOLOGY STUDY (M7)
% ========================================

clc;
clear;
close all;
addpath(genpath(pwd));
rng(1);

fprintf('\n========================================\n');
fprintf('   MAIN TEST 7 — TOPOLOGY STUDY (M7)\n');
fprintf('========================================\n');

%% =========================
% PART 1 — SINGLE CASE
% =========================
fprintf('\n--- PART 1: SINGLE DISCONNECTED CASE ---\n');

res_single = run_topology_study();

% --- GEOMETRY CLASSIFICATION ---
topology_flag = 1; % 🔥 CRITICAL (M7)

class_single = classify_geometry(0, 0, topology_flag);

fprintf('\nSUMMARY (SINGLE CASE)\n');
fprintf('Pf_MC   = %.6e\n', res_single.Pf_MC);
fprintf('Pf_FORM = %.6e\n', res_single.Pf_FORM);
fprintf('Pf_SORM = %.6e\n', res_single.Pf_SORM);

fprintf('Error FORM = %.4f\n', res_single.error_FORM);
fprintf('Error SORM = %.4f\n', res_single.error_SORM);

fprintf('FORM region = %s\n', res_single.region_selected);
fprintf('Geometry class = %s\n', class_single);

plot_disconnected_domain(res_single.problem, res_single.U_star);

fprintf('\n=== VISUALIZATION (M7 — SINGLE CASE) ===\n');
plot_limit_state_2d(res_single.problem, res_single.U_star);

%% =========================
% PART 2 — TOPOLOGY TRANSITION
% =========================
fprintf('\n--- PART 2: TOPOLOGY TRANSITION (c variation) ---\n');

results_c = run_topology_c_study();

n = length(results_c);

c_vals = zeros(n,1);
Pf_mc = zeros(n,1);
Pf_form = zeros(n,1);
Pf_sorm = zeros(n,1);
err_form = zeros(n,1);
err_sorm = zeros(n,1);

geom_class = strings(n,1); 

for i = 1:n
    c_vals(i) = results_c(i).c;

    Pf_mc(i) = results_c(i).Pf_MC;
    Pf_form(i) = results_c(i).Pf_FORM;
    Pf_sorm(i) = results_c(i).Pf_SORM;

    err_form(i) = results_c(i).error_FORM;
    err_sorm(i) = results_c(i).error_SORM;

    % --- CLASSIFICATION ---
    geom_class(i) = classify_geometry(0, 0, 1); % always disconnected
end

%% =========================
% SORT
% =========================
[c_vals, idx] = sort(c_vals);

Pf_mc = Pf_mc(idx);
Pf_form = Pf_form(idx);
Pf_sorm = Pf_sorm(idx);

err_form = err_form(idx);
err_sorm = err_sorm(idx);

geom_class = geom_class(idx); 

if i == 1 || i == round(n/2) || i == n
    fprintf('Plotting topology case c = %.2f\n', results_c(i).c);
    plot_limit_state_2d(results_c(i).problem, results_c(i).U_star);
end

%% =========================
% PRINT TABLE
% =========================
fprintf('\n========================================\n');
fprintf('   TOPOLOGY TRANSITION TABLE\n');
fprintf('========================================\n');

fprintf(' c     Pf_MC      Pf_FORM    Pf_SORM    err_FORM   err_SORM   class\n');

for i = 1:n
    fprintf('%4.2f  %10.3e  %10.3e  %10.3e  %8.3f   %8.3f   %s\n', ...
        c_vals(i), Pf_mc(i), Pf_form(i), Pf_sorm(i), ...
        err_form(i), err_sorm(i), ...
        geom_class(i));
end

%% =========================
% PLOT 1 — ERROR vs c
% =========================
figure;
plot(c_vals, err_form, 'o-', 'LineWidth', 2); hold on;
plot(c_vals, err_sorm, 's-', 'LineWidth', 2);

xlabel('Separation c');
ylabel('Relative Error');
title('Topology Effect: Error vs Separation');
legend('FORM','SORM','Location','best');
grid on;

%% =========================
% PLOT 2 — Pf comparison
% =========================
figure;
semilogy(c_vals, Pf_mc, 'k-o', 'LineWidth', 2); hold on;
semilogy(c_vals, Pf_form, 'r-s', 'LineWidth', 2);
semilogy(c_vals, Pf_sorm, 'b-d', 'LineWidth', 2);

xlabel('Separation c');
ylabel('Failure Probability (log scale)');
title('Pf Comparison vs Separation');
legend('MC','FORM','SORM','Location','best');
grid on;

%% =========================
% PLOT 3 — MISSED PROBABILITY
% =========================
figure;
plot(c_vals, 1 - Pf_form./Pf_mc, 'o-', 'LineWidth', 2); hold on;
plot(c_vals, 1 - Pf_sorm./Pf_mc, 's-', 'LineWidth', 2);

xlabel('Separation c');
ylabel('Missed Probability Fraction');
title('How Much Probability FORM/SORM Miss');
legend('FORM','SORM');
grid on;

%% =========================
% PLOT 4 — CONTRIBUTION SPLIT
% =========================
Pf_left  = [results_c.Pf_left];
Pf_right = [results_c.Pf_right];

figure;
plot(c_vals, Pf_left, 'o-', 'LineWidth', 2); hold on;
plot(c_vals, Pf_right, 's-', 'LineWidth', 2);

xlabel('Separation c');
ylabel('Contribution to Pf');
title('Failure Probability Split Between Regions');
legend('Left region','Right region');
grid on;

%% =========================
% PLOT 5 — CAPTURE RATIO
% =========================
capture = [results_c.capture_ratio];

figure;
plot(c_vals, capture, 'o-', 'LineWidth', 2);

xlabel('Separation c');
ylabel('Captured Fraction of Dominant Region');
title('FORM Accuracy on Single Failure Region');
grid on;

%% =========================
% INTERPRETATION
% =========================
fprintf('\n========================================\n');
fprintf('   TOPOLOGY INTERPRETATION\n');
fprintf('========================================\n');

for i = 1:n

    if err_form(i) < 0.1
        regime = 'CONNECTED (FORM OK)';
    elseif err_form(i) < 0.5
        regime = 'TRANSITION';
    else
        regime = 'DISCONNECTED (FORM FAILS)';
    end

    fprintf('c=%.2f | err_FORM=%.3f | class=%s | %s\n', ...
        c_vals(i), err_form(i), geom_class(i), regime);
end

fprintf('\n=== M7 FULL STUDY COMPLETE ===\n');