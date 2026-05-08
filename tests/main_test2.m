clc;
clear;

addpath(genpath(pwd));
rng(1);

fprintf('============================\n');
fprintf('        M2 FORM TEST\n');
fprintf('============================\n\n');

%% ============================
% TEST 1 — CANTILEVER
%% ============================

fprintf('--- CANTILEVER ---\n');

problem = get_problem_cantilever();
result = run_form(problem);

fprintf('\nRESULT:\n');
fprintf('Converged: %d\n', result.converged);
fprintf('Beta: %.6f\n', result.beta);
fprintf('Pf: %.6e\n', result.Pf);
fprintf('Evaluations: %d\n', result.neval);
fprintf('Runtime: %.4f s\n\n', result.runtime);

%% ============================
% TEST 2 — SYNTHETIC
%% ============================

fprintf('--- SYNTHETIC ---\n');

problem = get_problem_synthetic();
result = run_form(problem);

fprintf('\nRESULT:\n');
fprintf('Converged: %d\n', result.converged);
fprintf('Beta: %.6f\n', result.beta);
fprintf('Pf: %.6e\n', result.Pf);
fprintf('Evaluations: %d\n', result.neval);
fprintf('Runtime: %.4f s\n\n', result.runtime);

problem = get_problem_cantilever();
result = run_form(problem);

fprintf('\n--- M2 SANITY CHECK ---\n');

% 1. Design point consistency
g_star = evaluate_limit_state_u(problem, result.U_star);
fprintf('g(U*): %.3e\n', g_star);

% 2. Reliability identity
beta_check = norm(result.U_star);
fprintf('||U*||: %.6f vs beta: %.6f\n', beta_check, result.beta);

% 3. Alpha normalization
alpha_norm = norm(result.alpha);
fprintf('||alpha||: %.6f\n', alpha_norm);

% 4. Gradient direction sanity
dot_val = abs(dot(result.grad_g, result.alpha)) / norm(result.grad_g);
fprintf('dot(grad, alpha): %.6f\n', dot_val);

% 5. Iteration behavior (basic check)
fprintf('Iterations: %d\n', length(result.history.g));