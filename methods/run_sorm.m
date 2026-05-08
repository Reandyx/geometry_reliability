function result = run_sorm(problem)

    tic;
    reset_eval_count();

    res_form = run_form(problem);

    beta   = res_form.beta;
    U_star = res_form.U_star;
    grad_g = res_form.grad_g;
    Pf_form = res_form.Pf;

    result = create_result_struct();

    result.method     = 'SORM';
    result.beta       = beta;
    result.U_star     = U_star;
    result.neval      = eval_counter('get', 0);
    result.runtime    = toc;
    result.converged  = res_form.converged;
    result.u_star = result.U_star;

    if ~isfield(result, 'kappa')
        result.kappa = NaN;
    end
    
    tol = 1e-8;
    active_dims = sum(abs(grad_g) > tol);

    if active_dims <= 1
        warning('SORM skipped: degenerate (1D) problem detected. Using FORM.');

        % In 1D → FORM is exact → SORM = FORM
        result.Pf    = Pf_form;
        result.kappa = 0;

        return;
    end

    H = compute_hessian_u(problem, U_star);

    kappa = compute_principal_curvatures(H, grad_g);

    % ===== DEBUG BLOCK =====
    alpha = grad_g / norm(grad_g);

    % Tangent (2D case)
    t = [-alpha(2); alpha(1)];

    % Store curvature
    result.kappa = kappa;

    % Breitung correction 
    prod_term = 1;

    for i = 1:length(kappa)

        val = 1 + beta * kappa(i);

        % --- SORM breakdown handling ---
        if val <= 0
            warning('SORM breakdown: negative curvature term');

            result.Pf = NaN;
            return;
        end

        prod_term = prod_term * val^(-0.5);
    end

    Pf_sorm = Pf_form * prod_term;

    % --- Final assignment ---
    result.Pf = Pf_sorm;

    % --- Safety checks ---
    assert(isfield(result, 'Pf'));
    assert(isfield(result, 'beta'));
    assert(isfield(result, 'U_star'));

end