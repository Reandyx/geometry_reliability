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
        %fprintf('Iter %d: g = %.3e, ||U|| = %.4f\n', k, g, norm(U));
        
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
        % ==============================
        % SOFT-MIN DIAGNOSTIC
        % ==============================
        % if isfield(problem, 'debug_softmin') && problem.debug_softmin
        % 
        %     U_star = U;
        % 
        %     b = problem.debug_params.b;
        %     c = problem.debug_params.c;
        %     k_soft = problem.debug_params.k;
        % 
        %     g1 = b - sqrt(U_star(1)^2 + U_star(2)^2);
        %     g2 = b - sqrt((U_star(1)-c)^2 + U_star(2)^2);
        % 
        %     w1 = exp(-k_soft*g1);
        %     w2 = exp(-k_soft*g2);
        %     w1 = w1 / (w1 + w2);
        % 
        %     fprintf('[DEBUG softmin] k=%.1f | w1=%.3f\n', k_soft, w1);
        % 
        % end
    end
    
    % --- Final design point ---
    U_star = U;
    
    % Reliability index
    beta = norm(U_star);
    
    % Failure probability
    Pf = normcdf(-beta);
    
    % --- Final gradient ---
    grad_g_final = compute_gradient_u(problem, U_star);

    % --- Consistency check ---
    g_star = evaluate_limit_state_u(problem, U_star);
    g_star = g_star(1);   
    
    if abs(g_star) > 1e-5
        warning('FORM: Design point not exactly on limit-state surface');
    end

% ===== FULL FORM DEBUG =====
% fprintf('\n--- FORM FINAL DEBUG ---\n');
% 
% fprintf('U* = [%.6f, %.6f]\n', U_star(1), U_star(2));
% 
% fprintf('grad = [%.6e, %.6e]\n', grad_g_final(1), grad_g_final(2));
% fprintf('||grad|| = %.6e\n', norm(grad_g_final));
% 
% % alignment check (CRITICAL)
% cos_angle = dot(grad_g_final, U_star) / (norm(grad_g_final)*norm(U_star) + 1e-12);
% fprintf('cos(angle grad, U) = %.6f\n', cos_angle);
% 
% % limit-state check
% fprintf('g(U*) = %.6e\n\n', g_star);

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