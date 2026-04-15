clc;
clear;
close all;
addpath(genpath(pwd));
rng(1);

fprintf('========================================\n');
fprintf('   MAIN TEST 8.4 — IMPORTANCE SAMPLING\n');
fprintf('========================================\n\n');

% ========================================
% SETTINGS
% ========================================

N = 1e6;
c_vals = [0.1, 0.5, 1.0, 2.0, 3.0, 4.0];
n_cases = length(c_vals);

% Storage
beta_vals = zeros(n_cases,1);

Pf_MC   = zeros(n_cases,1);
Pf_IS   = zeros(n_cases,1);
Pf_FORM = zeros(n_cases,1);

CoV_MC  = zeros(n_cases,1);
CoV_IS  = zeros(n_cases,1);

fail_MC = zeros(n_cases,1);
fail_IS = zeros(n_cases,1);

gain    = zeros(n_cases,1);
ESS_rat = zeros(n_cases,1);

fprintf('--- RUNNING CASES ---\n\n');

% ========================================
% LOOP OVER CASES
% ========================================

for i = 1:n_cases

    c = c_vals(i);
    fprintf('--- CASE c = %.2f ---\n', c);

    % --- Problem ---
    problem = get_problem_synthetic_rare(3.0, c);

    % ========================================
    % FORM
    % ========================================
    res_form = run_form(problem);

    Pf_form = res_form.Pf;
    beta    = res_form.beta;

    % --- SAFE DESIGN POINT EXTRACTION ---
    if isfield(res_form, 'U_star')
        u_star = res_form.U_star;
    elseif isfield(res_form, 'u_star')
        u_star = res_form.u_star;
    else
        error('FORM does not return design point');
    end

    % ========================================
    % FORCE CORRECT GEOMETRY
    % ========================================
    if norm(u_star) > 0
        u_star = beta * (u_star / norm(u_star));
    else
        warning('u_star is zero — fallback used');
        u_star = beta * ones(size(u_star)) / sqrt(length(u_star));
    end

    % ========================================
    % MONTE CARLO
    % ========================================
    res_mc = run_mcs(problem, N);

    Pf_mc   = res_mc.Pf;

    % Safe CoV extraction
    if isfield(res_mc, 'CoV')
        CoV_mc = res_mc.CoV;
    else
        CoV_mc = res_mc.history.CoV;
    end

    fail_mc = res_mc.history.failures;

    % ========================================
    % IMPORTANCE SAMPLING
    % ========================================
    [Pf_is, CoV_is, res_is] = run_is(problem, N, u_star);

    % ========================================
    % STORE RESULTS
    % ========================================
    beta_vals(i) = beta;

    Pf_MC(i)   = Pf_mc;
    Pf_IS(i)   = Pf_is;
    Pf_FORM(i) = Pf_form;

    CoV_MC(i)  = CoV_mc;
    CoV_IS(i)  = CoV_is;

    fail_MC(i) = fail_mc;
    fail_IS(i) = res_is.failures;

    ESS_rat(i) = res_is.ESS_ratio;

    % Efficiency gain
    if CoV_is > 0 && isfinite(CoV_mc)
        gain(i) = (CoV_mc / CoV_is)^2;
    else
        gain(i) = Inf;
    end

    % ========================================
    % PRINT
    % ========================================
    fprintf('Pf_MC   = %.3e | CoV = %.3f\n', Pf_mc, CoV_mc);
    fprintf('Pf_IS   = %.3e | CoV = %.3f\n', Pf_is, CoV_is);
    fprintf('Pf_FORM = %.3e\n', Pf_form);
    fprintf('ESS = %.2f%% | Gain = %.2f\n\n', 100*ESS_rat(i), gain(i));

end

% ========================================
% SUMMARY TABLE
% ========================================

fprintf('========================================\n');
fprintf('   SUMMARY TABLE\n');
fprintf('========================================\n');

fprintf(' beta   Pf_MC      Pf_IS      Pf_FORM    CoV_MC   CoV_IS   gain   ESS   fail_MC  fail_IS\n');

for i = 1:n_cases
    fprintf('%5.2f  %9.3e  %9.3e  %9.3e   %6.3f   %6.3f  %6.2f  %5.2f   %6d   %6d\n', ...
        beta_vals(i), Pf_MC(i), Pf_IS(i), Pf_FORM(i), ...
        CoV_MC(i), CoV_IS(i), gain(i), ESS_rat(i), ...
        fail_MC(i), fail_IS(i));
end

% ========================================
% SORT BY BETA 
% ========================================
[beta_vals, idx] = sort(beta_vals);

Pf_MC   = Pf_MC(idx);
Pf_IS   = Pf_IS(idx);
Pf_FORM = Pf_FORM(idx);

CoV_MC  = CoV_MC(idx);
CoV_IS  = CoV_IS(idx);

gain    = gain(idx);
ESS_rat = ESS_rat(idx);

% ========================================
% FIGURE 1 — Pf comparison
% ========================================
figure;
semilogy(beta_vals, Pf_MC, 'ko-', 'LineWidth', 2, 'MarkerSize', 6); hold on;
semilogy(beta_vals, Pf_IS, 'bs-', 'LineWidth', 2, 'MarkerSize', 6);
semilogy(beta_vals, Pf_FORM, 'rd-', 'LineWidth', 2, 'MarkerSize', 6);

grid on;
xlabel('\beta');
ylabel('Failure Probability (log)');
title('Failure Probability Comparison');

legend('MC','IS','FORM','Location','southwest');

% ========================================
% FIGURE 2 — CoV comparison
% ========================================
figure;
plot(beta_vals, CoV_MC, 'ko-', 'LineWidth', 2, 'MarkerSize', 6); hold on;
plot(beta_vals, CoV_IS, 'bs-', 'LineWidth', 2, 'MarkerSize', 6);

grid on;
xlabel('\beta');
ylabel('Coefficient of Variation');
title('Variance Comparison');

legend('MC','IS');

% ========================================
% FIGURE 3 — Efficiency Gain
% ========================================
figure;
plot(beta_vals, gain, 'bd-', 'LineWidth', 2, 'MarkerSize', 6);

grid on;
xlabel('\beta');
ylabel('Gain');
title('Efficiency Gain of Importance Sampling');

% ========================================
% FIGURE 4 — ESS
% ========================================
figure;
plot(beta_vals, ESS_rat, 'ms-', 'LineWidth', 2, 'MarkerSize', 6);

grid on;
xlabel('\beta');
ylabel('ESS / N');
title('Effective Sample Size Ratio');

% ========================================
% INTERPRETATION
% ========================================

fprintf('\n========================================\n');
fprintf('   QUICK INTERPRETATION\n');
fprintf('========================================\n');

for i = 1:n_cases

    if gain(i) > 10
        msg = 'Strong variance reduction';
    elseif gain(i) > 2
        msg = 'Moderate improvement';
    else
        msg = 'Weak IS benefit';
    end

    fprintf('beta=%.2f | gain=%.2f | ESS=%.3f | %s\n', ...
        beta_vals(i), gain(i), ESS_rat(i), msg);
end

fprintf('\n=== M8.4 COMPLETE ===\n');