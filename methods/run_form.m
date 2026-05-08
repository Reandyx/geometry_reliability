function result = run_form(problem)
% RUN_FORM
% Stabilized First Order Reliability Method (HL-RF)
% Includes damping, line search, and strict validity checks

tic;
reset_eval_count();

% Initialization
dim = problem.dimension;

if isfield(problem, 'U0')
    U = problem.U0;
else
    U = ones(dim,1) * 0.5;
end

max_iter = 100;
tol_u = 1e-6;
tol_g = 1e-6;  %tol_g = 1e-3; can be used if data to noisy in real benchmarks

converged = false;

% History
history.U = [];
history.g = [];
history.normU = [];
history.step = [];
history.lambda = [];

for k = 1:max_iter
    
    % --- Evaluate limit-state ---
    g = evaluate_limit_state_u(problem, U);
    g = g(1);
    
    % --- Gradient ---
    grad_g = compute_gradient_u(problem, U);
    norm_grad = norm(grad_g);
    
    % Safeguard
    if norm_grad < 1e-12
        warning('FORM: Zero gradient encountered');
        break;
    end
    
    % HL-RF update
    U_hlrf = ((grad_g' * U - g) / (norm_grad^2)) * grad_g;
    
    % --- Step control ---
    step = U_hlrf - U;
    if norm(step) > 5
        step = 5 * step / norm(step);
    end
    U_hlrf = U + step;
    
    % LINE SEARCH (MERIT FUNCTION)
    lambda = 1.0;
    merit = 0.5 * norm(U)^2 + abs(g);
    
    while lambda > 1e-4
    
        U_trial = U + lambda * (U_hlrf - U);
        
        g_trial = evaluate_limit_state_u(problem, U_trial);
        g_trial = g_trial(1);
        
        merit_trial = 0.5 * norm(U_trial)^2 + 10 * abs(g_trial);
        
        if merit_trial < merit
            break;
        end
        
        lambda = lambda / 2;
    end
    
    % If line search failed → reject
    if lambda <= 1e-4
        warning('FORM: Line search failed — accepting step');
        U_new = U_hlrf;   % fallback
    else
        U_new = U_trial;
    end
    
    % Store history
    history.U(k,:) = U';
    history.g(k) = g;
    history.normU(k) = norm(U);
    history.step(k) = norm(U_new - U);
    history.lambda(k) = lambda;

    % Convergence check (FIXED)
    g_new = evaluate_limit_state_u(problem, U_new);
    g_new = g_new(1);
    
    if norm(U_new - U) < tol_u && abs(g_new) < tol_g
        converged = true;
        U = U_new;
        break;
    end
    
    % Divergence safeguard
    if norm(U_new) > 1e6
        warning('FORM: Divergence detected');
        break;
    end
    
    U = U_new;
end

n_iterations = k; 

if ~converged
    warning('FORM: Did not converge within max_iter');
end

% Final design point
U_star = U;

beta = norm(U_star);
if beta < 1e-12
    Pf = 0.5;
else
    Pf = normcdf(-beta);
end

% --- Final gradient ---
grad_g_final = compute_gradient_u(problem, U_star);
H = compute_hessian_u(problem, U_star);     %FORM conditioning evaluation AT at converged MPP 
cond_H = cond(H);

if ~isfinite(cond_H)
    cond_H = Inf;
end

% --- Residual ---
g_star = evaluate_limit_state_u(problem, U_star);
g_star = g_star(1);

% Validity check
if abs(g_star) > 1e-3
    warning('FORM: weak convergence (accepted)');
end

% Result struct
result = create_result_struct();

result.method = 'FORM';
result.beta = beta;
result.Pf = Pf;
result.U_star = U_star;
result.alpha = grad_g_final / norm(grad_g_final);
result.grad_g = grad_g_final;
result.neval = eval_counter('get', 0);
result.runtime = toc;
result.converged = converged;
result.history = history;
result.residual = abs(g_star);
result.u_star = result.U_star;
result.kappa = NaN;

if converged
    result.flag = 0;
else
    result.flag = 1;
end

result.n_iterations = n_iterations;
result.H = H;
result.cond_H = cond_H;
result.grad_norm = norm(grad_g_final);

% Safety checks
assert(isfield(result, 'Pf'));
assert(isfield(result, 'beta'));
assert(isfield(result, 'U_star'));

end
