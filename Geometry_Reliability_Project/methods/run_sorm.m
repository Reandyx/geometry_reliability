function result = run_sorm(problem)

    tic;
    reset_eval_count();

    % --- Step 1: FORM ---
    res_form = run_form(problem);

    beta   = res_form.beta;
    U_star = res_form.U_star;
    grad_g = res_form.grad_g;
    Pf_form = res_form.Pf;

    % --- Step 2: Initialize result ---
    result = create_result_struct();

    result.method     = 'SORM';
    result.beta       = beta;
    result.U_star     = U_star;
    result.neval      = eval_counter('get', 0);
    result.runtime    = toc;
    result.converged  = res_form.converged;

    % ============================================================
    % --- Step 3: Degeneracy detection 
    % ============================================================
    tol = 1e-8;
    active_dims = sum(abs(grad_g) > tol);

    if active_dims <= 1
        warning('SORM skipped: degenerate (1D) problem detected. Using FORM.');

        % In 1D → FORM is exact → SORM = FORM
        result.Pf    = Pf_form;
        result.kappa = 0;

        return;
    end

    % ============================================================
    % --- Step 4: Hessian ---
    % ============================================================
    H = compute_hessian_u(problem, U_star);

    % ============================================================
    % --- Step 5: Principal curvatures ---
    % ============================================================
    kappa = compute_principal_curvatures(H, grad_g);

    % ===== DEBUG BLOCK =====
    alpha = grad_g / norm(grad_g);

    % Tangent (2D case)
    t = [-alpha(2); alpha(1)];

    fprintf('\n--- SORM GEOMETRY DEBUG ---\n');
    disp('U* ='); disp(U_star);
    disp('grad_g ='); disp(grad_g);
    disp('H ='); disp(H);

    disp('alpha (normal) ='); disp(alpha);
    disp('tangent t ='); disp(t);

    fprintf('t^T H t = %.6e\n', t' * H * t);
    fprintf('||grad|| = %.6e\n', norm(grad_g));

    if ~isempty(kappa)
        disp('kappa (principal) ='); disp(kappa);
    else
        disp('kappa is EMPTY');
    end

    % Store curvature
    result.kappa = kappa;

    % ============================================================
    % --- Step 6: Breitung correction ---
    % ============================================================
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

    % ============================================================
    % --- DEBUG OUTPUT ---
    % ============================================================
    fprintf('\n--- DEBUG SORM ---\n');
    disp('U* ='); disp(U_star);
    disp('grad_g ='); disp(grad_g);
    disp('Hessian ='); disp(H);
    disp('Principal curvatures kappa ='); disp(kappa);

    fprintf('beta = %.6f\n', beta);

    for i = 1:length(kappa)
        fprintf('1 + beta*kappa(%d) = %.6f\n', i, 1 + beta * kappa(i));
    end

    % --- Safety checks ---
    assert(isfield(result, 'Pf'));
    assert(isfield(result, 'beta'));
    assert(isfield(result, 'U_star'));

end