clc;
clear;
addpath(genpath(pwd));
rng(1);

disp('=== M1 SYSTEM TEST START ===');

%% =========================
% 1. LOAD PROBLEM
%% =========================
problem = get_problem_cantilever();

nSamples = 10000;

%% =========================
% 2. BASIC MCS TEST
%% =========================
disp('--- Monte Carlo Test ---');

X = sample_inputs(problem, nSamples);

tic;
g = evaluate_limit_state(problem, X);
Pf = mean(g <= 0);
runtime = toc;

beta = compute_beta_from_pf(Pf);

disp(['Pf:   ', num2str(Pf)]);
disp(['beta: ', num2str(beta)]);
disp(['runtime [s]: ', num2str(runtime)]);

%% =========================
% 3. TRANSFORMATION CONSISTENCY TEST
%% =========================
disp('--- Transformation Consistency Test ---');

[nX, d] = size(X);

X_back = zeros(nX, d);

for k = 1:nX
    xk = X(k, :)';          % (d x 1)

    uk = x_to_u(xk, problem);
    xk_back = u_to_x(uk, problem);

    X_back(k, :) = xk_back';   % back to row
end

err_transform = max(abs(X(:) - X_back(:)));

disp(['Transform error: ', num2str(err_transform)]);

%% =========================
% 4. U-SPACE EVALUATION TEST
%% =========================
disp('--- U-Space Evaluation Test ---');

g_u = zeros(nX, 1);

for k = 1:nX
    xk = X(k, :)';          % (d x 1)

    uk = x_to_u(xk, problem);
    g_u(k) = evaluate_limit_state_u(problem, uk);
end

err_u = max(abs(g(:) - g_u(:)));

disp(['U-space consistency error: ', num2str(err_u)]);

%% =========================
% 5. SYNTHETIC PROBLEM TEST
%% =========================
disp('--- Synthetic Problem Test ---');

problem_syn = get_problem_synthetic();

% baseline curvature
X_syn = sample_inputs(problem_syn, nSamples);
g_syn = evaluate_limit_state(problem_syn, X_syn);
Pf_syn1 = mean(g_syn <= 0);

disp(['Pf (a=1.0): ', num2str(Pf_syn1)]);

% modified curvature
problem_syn.a = 0.5;

X_syn = sample_inputs(problem_syn, nSamples);
g_syn = evaluate_limit_state(problem_syn, X_syn);
Pf_syn2 = mean(g_syn <= 0);

disp(['Pf (a=0.5): ', num2str(Pf_syn2)]);

%% =========================
% 6. RESULT STRUCT TEST
%% =========================
disp('--- Result Struct Test ---');

result = create_result_struct();
result.method = 'MCS';
result.Pf = Pf;
result.beta = beta;
result.neval = nSamples;
result.runtime = runtime;
result.g_eval = g;
result.converged = true;

disp(result);

%% =========================
% 7. FINAL STATUS
%% =========================
disp('=== M1 SYSTEM TEST COMPLETE ===');