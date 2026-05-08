clc;
clear;
close all;
addpath(genpath(pwd));
rng(1);


fprintf('========================================\n');
fprintf(' VERIFICATION SUITE — SECTIONS 3–6\n');
fprintf('========================================\n\n');

results = struct();

%% ============================================================
% TEST 1 — GEOMETRY CLASSIFICATION (CURVATURE VS ERROR)
% ============================================================
fprintf('--- TEST 1: GEOMETRY CLASSIFICATION ---\n');

a_values = linspace(0.2, 5, 10);   % controls curvature
b = 4.0;                           % controls Pf level
N_mcs = 1e5;

geom_results = struct();

for i = 1:length(a_values)

    a = a_values(i);

    problem = get_problem_synthetic(a, b);

    % --- METHODS ---
    res_form = run_form(problem);
    res_mcs  = run_mcs(problem, N_mcs);

    % --- ERROR ---
    Pf_form = res_form.Pf;
    Pf_mcs  = res_mcs.Pf;

    if Pf_mcs == 0
        rel_error = NaN;
    else
        rel_error = abs(Pf_form - Pf_mcs) / Pf_mcs;
    end

    % --- CURVATURE ---
    kappa = compute_curvature_at_design_point(problem, res_form.U_star);
    kappa_max = max(abs(kappa));

    % --- GEOMETRY INDEX ---
    gamma = res_form.beta * kappa_max;

    % --- STORE ---
    geom_results(i).a = a;
    geom_results(i).Pf_form = Pf_form;
    geom_results(i).Pf_mcs = Pf_mcs;
    geom_results(i).rel_error = rel_error;
    geom_results(i).beta = res_form.beta;
    geom_results(i).kappa = kappa;
    geom_results(i).gamma = gamma;

    fprintf('a=%.2f | Pf_FORM=%.2e | Pf_MCS=%.2e | err=%.3f | gamma=%.3f\n', ...
        a, Pf_form, Pf_mcs, rel_error, gamma);
end

results.geometry = geom_results;


%% ============================================================
% TEST 2 — CURVATURE SWEEP (FORM vs SORM vs MCS)
% ============================================================
fprintf('\n--- TEST 2: CURVATURE SWEEP ---\n');

curv_results = struct();

for i = 1:length(a_values)

    a = a_values(i);

    problem = get_problem_synthetic(a, b);

    res_form = run_form(problem);
    res_sorm = run_sorm(problem);
    res_mcs  = run_mcs(problem, N_mcs);

    curv_results(i).a = a;
    curv_results(i).Pf_form = res_form.Pf;
    curv_results(i).Pf_sorm = res_sorm.Pf;
    curv_results(i).Pf_mcs = res_mcs.Pf;
    curv_results(i).beta = res_form.beta;
    curv_results(i).kappa = res_sorm.kappa;

    fprintf('a=%.2f | FORM=%.2e | SORM=%.2e | MCS=%.2e\n', ...
        a, res_form.Pf, res_sorm.Pf, res_mcs.Pf);
end

results.curvature = curv_results;


%% ============================================================
% TEST 3 — TOPOLOGY (DISCONNECTED FAILURE DOMAIN)
% ============================================================
fprintf('\n--- TEST 3: TOPOLOGY (DISCONNECTED) ---\n');

problem = get_problem_disconnected();

res_form = run_form(problem);
res_mcs  = run_mcs(problem, N_mcs);

Pf_form = res_form.Pf;
Pf_mcs  = res_mcs.Pf;

if Pf_mcs == 0
    rel_error = NaN;
else
    rel_error = abs(Pf_form - Pf_mcs) / Pf_mcs;
end

fprintf('Disconnected case:\n');
fprintf('FORM Pf = %.3e\n', Pf_form);
fprintf('MCS  Pf = %.3e\n', Pf_mcs);
fprintf('Relative error = %.3f\n', rel_error);

topo_results.Pf_form = Pf_form;
topo_results.Pf_mcs  = Pf_mcs;
topo_results.rel_error = rel_error;

results.topology = topo_results;


%% ============================================================
% TEST 4 — RARE EVENT STUDY (COST SCALING)
% ============================================================
fprintf('\n--- TEST 4: RARE EVENT STUDY ---\n');

beta_targets = [3, 4, 5, 6];

rare_results = struct();

for i = 1:length(beta_targets)

    beta_target = beta_targets(i);

    % calibrate threshold b for target beta
    problem = get_problem_synthetic(1.0, beta_target^2);

    % FORM for design point
    res_form = run_form(problem);

    % MCS
    N_test = 1e5;
    res_mcs = run_mcs(problem, N_test);

    % IS
    [Pf_is, CoV_is, is_info] = run_is(problem, 1e4, res_form.U_star);

    rare_results(i).beta_target = beta_target;
    rare_results(i).Pf_form = res_form.Pf;
    rare_results(i).Pf_mcs = res_mcs.Pf;
    rare_results(i).Pf_is  = Pf_is;
    rare_results(i).CoV_is = CoV_is;
    rare_results(i).ESS    = is_info.ESS;

    fprintf('beta=%.1f | FORM=%.2e | MCS=%.2e | IS=%.2e | CoV_IS=%.2f\n', ...
        beta_target, res_form.Pf, res_mcs.Pf, Pf_is, CoV_is);
end

results.rare_event = rare_results;


%% ============================================================
% SAVE RESULTS
% ============================================================
save(fullfile('results', 'results_verification_suite_main_test9.mat'), 'results');

fprintf('\n========================================\n');
fprintf(' ALL TESTS COMPLETED\n');
fprintf('========================================\n');


%% ============================================================
% Helper function
% ============================================================
function kappa = compute_curvature_at_design_point(problem, u_star)

grad = compute_gradient_u(problem, u_star);
H    = compute_hessian_u(problem, u_star);

alpha = grad / norm(grad);
P = eye(length(u_star)) - alpha * alpha';

H_proj = P * H * P;

kappa = eig(H_proj) / norm(grad);

end