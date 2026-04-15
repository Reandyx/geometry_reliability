clc;
clear;
close all;
addpath(genpath(pwd));
rng(1);

fprintf('\n========================================\n');
fprintf(' CONTROLLED VALIDATION SUITE — TEST 10\n');
fprintf('========================================\n');

results = struct();

%% ============================================================
% BLOCK 4 — FIXED β, VARYING CURVATURE
% ============================================================
fprintf('\n--- BLOCK 4: FIXED BETA (CURVATURE) ---\n');

target_beta = 2.0;
Pf_target = normcdf(-target_beta);

a_values = linspace(0.01, 0.5, 20);
c = 0.0;

for i = 1:length(a_values)

    a = a_values(i);

    % --- calibrate ---
    b = calibrate_b(a, Pf_target, 'global');
    problem = get_problem_synthetic(a, b, c);

    problem.grad_g = @(U) compute_gradient_u(problem, U);
    problem.hess_g = @(U) compute_hessian_u(problem, U);

    % --- run methods ---
    res_form = run_form(problem);
    res_sorm = run_sorm(problem);
    res_mcs  = run_mcs(problem, 1e5);

    % ============================================================
    % CURVATURE (κ) FROM HESSIAN
    % ============================================================
    U_star = res_form.U_star(:)';

    grad = problem.grad_g(U_star);
    grad = grad(:);

    H = problem.hess_g(U_star);

    alpha = grad / norm(grad);

    P = eye(length(U_star)) - alpha * alpha';

    H_t = P * H * P;

    eigvals = eig(H_t);

    % remove near-zero eigenvalue (normal direction)
    kappa = eigvals(abs(eigvals) > 1e-6);

    % ============================================================
    % BREITUNG SORM (CONSISTENT)
    % ============================================================
    if isempty(kappa)
        Pf_breitung = res_form.Pf;
        gamma = 0;
    else
        breitung_factor = prod((1 + res_form.beta * kappa).^(-0.5));
        Pf_breitung = normcdf(-res_form.beta) * breitung_factor;

        gamma = res_form.beta * max(abs(kappa));
    end

    % ============================================================
    % ERRORS
    % ============================================================
    err_form = abs(res_form.Pf - res_mcs.Pf) / res_mcs.Pf;
    err_sorm = abs(res_sorm.Pf - res_mcs.Pf) / res_mcs.Pf;
    err_breitung = abs(Pf_breitung - res_mcs.Pf) / res_mcs.Pf;

    % ============================================================
    % STORE
    % ============================================================
    results.curvature(i).a = a;
    results.curvature(i).beta = res_form.beta;
    results.curvature(i).gamma = gamma;
    results.curvature(i).Pf_form = res_form.Pf;
    results.curvature(i).Pf_sorm = res_sorm.Pf;
    results.curvature(i).Pf_mcs = res_mcs.Pf;
    results.curvature(i).err_form = err_form;
    results.curvature(i).err_sorm = err_sorm;
    results.curvature(i).Pf_breitung = Pf_breitung;
    results.curvature(i).err_breitung = err_breitung;
    results.curvature(i).kappa = kappa;

    % ============================================================
    % LOG
    % ============================================================
    fprintf('a=%.2f | beta=%.2f | gamma=%.3f | err_FORM=%.3f | err_SORM=%.3f | err_Breit=%.3f\n', ...
        a, res_form.beta, gamma, err_form, err_sorm, err_breitung);
end

% gamma_vals = [results.curvature.gamma];
% disp(gamma_vals)

%% ============================================================
% CLASSIFICATION SUMMARY
% ============================================================
fprintf('\n--- CLASSIFICATION SUMMARY ---\n');

gamma_vals = [results.curvature.gamma];
err_vals   = [results.curvature.err_form];

[gamma_sorted, idx] = sort(gamma_vals);
err_sorted = err_vals(idx);

% --- transition detection ---
err_smooth = movmean(err_sorted, 3);
d_err = diff(err_smooth);

if ~isempty(d_err)
    [~, i_transition] = max(d_err);
    gamma_transition = gamma_sorted(i_transition);
    fprintf('\nEstimated curvature threshold: gamma ≈ %.3f\n', gamma_transition);
else
    gamma_transition = NaN;
end

% --- fixed classes ---
gamma_I_max = 0.2;
gamma_II_max = 0.6;

class_I   = err_vals(gamma_vals < gamma_I_max);
class_II  = err_vals(gamma_vals >= gamma_I_max & gamma_vals <= gamma_II_max);
class_III = err_vals(gamma_vals > gamma_II_max);

fprintf('Class I   (gamma < %.1f): mean err = %.3f\n', gamma_I_max, mean(class_I));
fprintf('Class II  (%.1f–%.1f):    mean err = %.3f\n', gamma_I_max, gamma_II_max, mean(class_II));
fprintf('Class III (gamma > %.1f): mean err = %.3f\n', gamma_II_max, mean(class_III));

% --- debug print ---
fprintf('\n--- RAW DATA (gamma vs err_FORM) ---\n');
for i = 1:length(gamma_sorted)
    fprintf('gamma = %.3f | err = %.3f\n', gamma_sorted(i), err_sorted(i));
end

%% ============================================================
% BLOCK 5A — CONTROLLED TOPOLOGY (DISCONNECTED FAILURE)
% ============================================================
fprintf('\n--- BLOCK 5A: TOPOLOGY (DISCONNECTED FAILURE) ---\n');

target_beta = 3.0;
Pf_target = normcdf(-target_beta);

c_values = linspace(0.5, 3.0, 10);

for i = 1:length(c_values)

    c = c_values(i);

    % --- calibrate ---
    b = calibrate_b_disconnected_c(c, Pf_target);
    problem = get_problem_disconnected_c(b, c);

    a = problem.a;

    % ============================================================
    % --- DEFINE BRANCH FUNCTIONS (SMOOTH + STABLE)
    % ============================================================
    g1_fun = @(U) (U(:,2).^2 + a*(U(:,1) - c).^2) - b;
    g2_fun = @(U) (U(:,2).^2 + a*(U(:,1) + c).^2) - b + 0.05;

    % --- FORM on branch 1 ---
    problem1 = problem;
    problem1.gfun = g1_fun;
    problem1.U0 = [c; 0];
    res1 = run_form(problem1);

    % --- FORM on branch 2 ---
    problem2 = problem;
    problem2.gfun = g2_fun;
    problem2.U0 = [-c; 0];
    res2 = run_form(problem2);

    % --- collect betas ---
    betas = [res1.beta; res2.beta];
    beta_form = min(betas);
    Pf_form = normcdf(-beta_form);

    % ============================================================
    % --- TRUE MODEL (SMOOTH MIN — IMPORTANT)
    % ============================================================
    k = 20;
    problem_full = problem;
    problem_full.gfun = @(U) ...
        -log(exp(-k*g1_fun(U)) + exp(-k*g2_fun(U))) / k;

    % --- MCS ---
    res_mcs = run_mcs(problem_full, 1e6);

    % ============================================================
    % --- ERROR
    % ============================================================
    err = abs(Pf_form - res_mcs.Pf) / res_mcs.Pf;

    % ============================================================
    % --- STORE
    % ============================================================
    results.topology(i).c = c;
    results.topology(i).beta = beta_form;
    results.topology(i).Pf_mcs = res_mcs.Pf;
    results.topology(i).Pf_form = Pf_form;
    results.topology(i).error = err;

    % ============================================================
    % --- LOGGING
    % ============================================================
    fprintf(['c=%.2f | beta=%.2f | Pf_FORM=%.2e | Pf_MCS=%.2e | err=%.2f\n'], ...
        c, beta_form, Pf_form, res_mcs.Pf, err);

end

fprintf('\n--- BLOCK 5B: SOFT TOPOLOGY SWEEP (k) ---\n');

target_beta = 3.0;
Pf_target = normcdf(-target_beta);

c_values = linspace(0.5, 3.0, 10);
k_values = [1, 5, 20, 50, 100];

for j = 1:length(k_values)

    k = k_values(j);
    fprintf('\n--- k = %.2f ---\n', k);

    for i = 1:length(c_values)

        c = c_values(i);

        % --- calibrate ---
        b = calibrate_b_softmin(c, Pf_target, k);

        beta1 = b;
        beta2 = b + 1.0;

        problem = get_problem_softmin(beta1, beta2, c, k);

        U0_list = [
            0, 0;
            3, 0;
            -3, 0;
            0, 3;
            0, -3
        ];

        betas = zeros(size(U0_list,1),1);

        for s = 1:size(U0_list,1)
            problem.U0 = U0_list(s,:)';
            res_tmp = run_form(problem);
            betas(s) = res_tmp.beta;
        end

        betas_sorted = sort(betas);

        beta_form = betas_sorted(1);
        Pf_form = normcdf(-beta_form);

        % --- UNIQUE MODES ---
        tol = 0.001;
        unique_betas = [];

        for s = 1:length(betas)
            if isempty(unique_betas) || all(abs(unique_betas - betas(s)) > tol)
                unique_betas(end+1) = betas(s);
            end
        end

        n_modes = length(unique_betas);
        beta_gap = max(unique_betas) - min(unique_betas);

        if beta_gap < 0.05
            topo_class = 1;
        elseif beta_gap < 0.5
            topo_class = 2;
        else
            topo_class = 3;
        end

        % --- MCS (HARD MIN for truth) ---
        g1_fun = @(U) beta1 - abs(U(:,1));
        g2_fun = @(U) beta2 - abs(U(:,2) + c);

        problem_mcs = problem;
        problem_calib.gfun = @(U) min(g1_fun(U), g2_fun(U));

        % IMPORTANT: MCS must use SAME soft function
        res_mcs = run_mcs(problem_mcs, 1e6);

        err = abs(Pf_form - res_mcs.Pf) / res_mcs.Pf;

        results.topology_soft(j).k = k;

        results.topology_soft(j).data(i).c = c;
        results.topology_soft(j).data(i).beta = beta_form;
        results.topology_soft(j).data(i).Pf_form = Pf_form;
        results.topology_soft(j).data(i).Pf_mcs = res_mcs.Pf;
        results.topology_soft(j).data(i).error = err;
        results.topology_soft(j).data(i).beta_gap = beta_gap;
        results.topology_soft(j).data(i).n_modes = n_modes;
        results.topology_soft(j).data(i).topo_class = topo_class;

        fprintf(['k=%.2f | c=%.2f | beta=%.2f | gap=%.3f | modes=%d | class=%d | err=%.3f\n'], ...
            k, c, beta_form, beta_gap, n_modes, topo_class, err);

    end
end

%% ============================================================
% BLOCK 3 — GEOMETRY CLASSIFICATION (ANCHORING)
% ============================================================
fprintf('\n--- BLOCK 3: GEOMETRY CLASSIFICATION (ANCHORING) ---\n');

target_beta = 3.0;
Pf_target = normcdf(-target_beta);

% --- single shape parameter sweep ---
a_values = linspace(0.05, 3.0, 15);

for i = 1:length(a_values)

    a = a_values(i);

    % ============================================================
    % --- CALIBRATE b (FIXED RELIABILITY LEVEL)
    % ============================================================
    b = calibrate_b(a, Pf_target, 'global');

    problem = get_problem_synthetic(a, b);

    % ============================================================
    % --- RUN METHODS
    % ============================================================
    res_form = run_form(problem);
    res_sorm = run_sorm(problem);
    res_mcs  = run_mcs(problem, 1e5);

    % ============================================================
    % --- CURVATURE (γ = β * max|κ|)
    % ============================================================
    if isempty(res_sorm.kappa)
        gamma = 0;
        kappa_max = 0;
    else
        kappa_max = max(abs(res_sorm.kappa));
        gamma = res_form.beta * kappa_max;
    end

    % ============================================================
    % --- ERROR (REFERENCE: MCS)
    % ============================================================
    err_form = abs(res_form.Pf - res_mcs.Pf) / res_mcs.Pf;

    % ============================================================
    % --- CLASSIFICATION (NUMERICALLY ANCHORED)
    % ============================================================
    if err_form < 0.05
        class = 1; % Class I
    elseif err_form < 0.30
        class = 2; % Class II
    else
        class = 3; % Class III
    end

    % ============================================================
    % --- STORE
    % ============================================================
    results.section3(i).a = a;
    results.section3(i).beta = res_form.beta;
    results.section3(i).gamma = gamma;
    results.section3(i).kappa = kappa_max;
    results.section3(i).err_form = err_form;
    results.section3(i).class = class;

    fprintf('a=%.2f | beta=%.2f | gamma=%.3f | kappa=%.3f | err=%.3f | class=%d\n', ...
        a, res_form.beta, gamma, kappa_max, err_form, class);

end

% ============================================================
% --- CLASS SUMMARY (KEY RESULT FOR PAPER)
% ============================================================
fprintf('\n--- CLASSIFICATION SUMMARY (NUMERICAL ANCHORING) ---\n');

gamma_vals = [results.section3.gamma];
err_vals   = [results.section3.err_form];

class_I   = err_vals(err_vals < 0.05);
class_II  = err_vals(err_vals >= 0.05 & err_vals < 0.30);
class_III = err_vals(err_vals >= 0.30);

fprintf('Class I   (near-linear):     err < 5%%    → mean = %.3f\n', mean(class_I));
fprintf('Class II  (curved):          5–30%%        → mean = %.3f\n', mean(class_II));
fprintf('Class III (high curvature):  >30%%         → mean = %.3f\n', mean(class_III));

% ============================================================
% --- TRANSITION DETECTION (γ threshold)
% ============================================================
[gamma_sorted, idx] = sort(gamma_vals);
err_sorted = err_vals(idx);

err_smooth = movmean(err_sorted, 3);
d_err = diff(err_smooth);

if ~isempty(d_err)
    [~, i_transition] = max(d_err);
    gamma_transition = gamma_sorted(i_transition);

    fprintf('\nEstimated transition threshold: gamma ≈ %.3f\n', gamma_transition);
end

% ============================================================
% --- RAW DATA (MANDATORY FOR PAPER)
% ============================================================
fprintf('\n--- RAW DATA (gamma vs err_FORM) ---\n');
for i = 1:length(gamma_sorted)
    fprintf('gamma = %.3f | err = %.3f\n', gamma_sorted(i), err_sorted(i));
end

%% ============================================================
%BLOCK 6 — RARE EVENT EFFICIENCY TABLE
%============================================================
fprintf('\n--- BLOCK 3: RARE EVENT TABLE ---\n');

beta_values = [3 4 5 6];
CoV_target = 0.10;

for i = 1:length(beta_values)

    beta_target = beta_values(i);

    problem = get_problem_synthetic_rare(beta_target);

    % --- FORM ---
    res_form = run_form(problem);
    Pf = res_form.Pf;

    % --- MCS required samples ---
    N_mcs = ceil(100 / Pf);   % CoV ≈ 0.1 rule

    % --- IS estimate ---
    N_is = 1e5;
    [Pf_IS, CoV_IS, ~] = run_is(problem, N_is, res_form.U_star);

    % --- required N for IS ---
    if CoV_IS > 0
        N_is_required = ceil(N_is * (CoV_IS / CoV_target)^2);
    else
        N_is_required = NaN;
    end

    % --- store ---
    results.rare(i).beta = beta_target;
    results.rare(i).Pf = Pf;
    results.rare(i).N_mcs = N_mcs;
    results.rare(i).N_is = N_is_required;
    results.rare(i).CoV_is = CoV_IS;

    fprintf('beta=%.1f | Pf=%.2e | N_MCS≈%.2e | N_IS≈%.2e | CoV_IS=%.2f\n', ...
        beta_target, Pf, N_mcs, N_is_required, CoV_IS);
end


%% ============================================================
% SAVE RESULTS
% ============================================================
if ~exist('results', 'dir')
    mkdir('results');
end

save(fullfile('results', 'results_main_test10.mat'), 'results');

fprintf('\n========================================\n');
fprintf(' TEST 10 COMPLETED\n');
fprintf('========================================\n');