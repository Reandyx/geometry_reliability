function result = run_form(problem)
% RUN_FORM
% First Order Reliability Method (HL-RF)

    tic;
    reset_eval_count();
    
    % Initialization
    dim = problem.dimension;
    U = zeros(dim,1);
    
    max_iter = 100;
    tol_u = 1e-6;
    tol_g = 1e-6;
    
    converged = false;
    
    % History
    history.U = [];
    history.g = [];
    history.normU = [];
    
    for k = 1:max_iter
        
        % --- Evaluate limit-state (force scalar) ---
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
        
        % --- Store history ---
        history.U(k,:) = U';
        history.g(k) = g;                 
        history.normU(k) = norm(U);
        
        % --- HL-RF update ---
        U_new = ((grad_g' * U - g) / (norm_grad^2)) * grad_g;
        
        % Debug
        fprintf('Iter %d: g = %.3e, ||U|| = %.4f\n', k, g, norm(U));
        
        % Convergence
        if norm(U_new - U) < tol_u && abs(g) < tol_g
            converged = true;
            U = U_new;
            break;
        end
        
        % Divergence
        if norm(U_new) > 1e6
            warning('FORM: Divergence detected');
            break;
        end
        
        U = U_new;
    end
    
    % --- Final design point ---
    U_star = U;
    
    % Reliability index
    beta = norm(U_star);
    
    % Failure probability
    Pf = normcdf(-beta);
    
    % --- Final gradient ---
    grad_g_final = compute_gradient_u(problem, U_star);
    
    % DEBUG CHECK
    fprintf('grad_g(2) norm = %.4e\n', abs(grad_g_final(2)));

    % --- Consistency check ---
    g_star = evaluate_limit_state_u(problem, U_star);
    g_star = g_star(1);   
    
    if abs(g_star) > 1e-5
        warning('FORM: Design point not exactly on limit-state surface');
    end
    
    % --- Result struct ---
    result = create_result_struct();
    
    result.method = 'FORM';
    result.beta = beta;
    result.Pf = Pf;
    result.U_star = U_star;
    result.alpha = U_star / beta;
    result.grad_g = grad_g_final;
    result.neval = eval_counter('get', 0);
    result.runtime = toc;
    result.converged = converged;
    result.history = history;
    
    % Safety checks
    assert(isfield(result, 'Pf'));
    assert(isfield(result, 'beta'));
    assert(isfield(result, 'U_star'));

end