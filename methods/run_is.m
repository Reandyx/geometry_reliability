function result = run_is(problem, u_star, N)

    tic;
    reset_eval_count();

    dim = problem.dimension;

    % --- IS distribution ---
    mu_is = u_star(:)';  

    sigma_scale = 1.5;
    Sigma_is = (sigma_scale^2) * eye(dim);

    U = mvnrnd(mu_is, Sigma_is, N);

    % ================================
    % TRANSFORM (CONSISTENT)
    % ================================
    X = u_to_x(U, problem);

    % ================================
    % LIMIT STATE (SINGLE SOURCE)
    % ================================
    g = evaluate_limit_state(problem, X);
    I = (g <= 0);

    % ================================
    % LOG WEIGHTS (STABLE)
    % ================================
    log_pdf_target = -0.5 * sum(U.^2, 2) - (dim/2)*log(2*pi);

    diff = U - mu_is;

    log_pdf_is = -0.5 * sum(diff.^2, 2) / (sigma_scale^2) ...
                 - (dim/2)*log(2*pi) ...
                 - dim*log(sigma_scale);

    log_w = log_pdf_target - log_pdf_is;

    % stabilization
    log_w = log_w - max(log_w);

    w = exp(log_w);

    % ================================
    % ESTIMATION
    % ================================
    w_sum = sum(w);

    if w_sum == 0
        warning('IS: all weights collapsed');
        Pf_IS = NaN;
        CoV_IS = NaN;
    else
        Pf_IS = sum(w .* I) / w_sum;

        var_IS = sum(w.^2 .* (I - Pf_IS).^2) / (w_sum^2);

        if Pf_IS < 1e-12
            CoV_IS = NaN;
        else
            CoV_IS = sqrt(var_IS) / Pf_IS;
        end
    end

    % ================================
    % DIAGNOSTICS
    % ================================
    den = sum(w.^2);
    
    if den == 0
        ESS = 0;
    else
        ESS = (sum(w)^2) / den;
    end
    
    ESS_ratio = ESS / N;

    fail_IS = sum(I);

    is_degenerate = (ESS_ratio < 0.05);
    is_no_failure = (fail_IS == 0);

    % ================================
    % PRINT
    % ================================
    fprintf('[IS] Pf = %.3e | CoV = %.3f | ESS = %.0f (%.1f%%) | Failures = %d\n', ...
        Pf_IS, CoV_IS, ESS, 100*ESS_ratio, fail_IS);

    if is_degenerate
        fprintf('[IS WARNING] Weight degeneracy detected (ESS < 5%%)\n');
    end

    if is_no_failure
        fprintf('[IS WARNING] No failure samples detected\n');
    end

    % ================================
    % RESULT STRUCT (R6 FIX)
    % ================================
    result = create_result_struct();

    result.method  = 'IS';
    result.Pf      = Pf_IS;

    if isnan(Pf_IS) || Pf_IS == 0
        result.beta = NaN;
    else
        result.beta = compute_beta_from_pf(Pf_IS);
    end

    result.neval   = N;
    result.runtime = toc;
    
    result.COV = CoV_IS;
    result.u_star  = u_star;
    result.kappa   = NaN;

    % diagnostics
    result.history.CoV        = CoV_IS;
    result.history.ESS        = ESS;
    result.history.ESS_ratio  = ESS_ratio;
    result.history.failures   = fail_IS;
    result.history.degenerate = is_degenerate;

end