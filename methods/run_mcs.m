function result = run_mcs(problem, N)
% RUN_MCS
% Crude Monte Carlo Simulation for reliability estimation

    tic;

    reset_eval_count();

    % --- Sampling ---
    X = sample_inputs(problem, N);

    % --- Limit state ---
    g_vals = evaluate_limit_state(problem, X);

    % --- Indicator ---
    indicator = (g_vals <= 0);

    % --- Probability ---
    Pf = mean(indicator);

    % --- Reliability index ---
    if Pf == 0
        beta = Inf;
    else
        beta = compute_beta_from_pf(Pf);
    end
    
    % --- Coefficient of Variation ---
    if Pf > 0
        CoV = sqrt((1 - Pf) / (N * Pf));
    else
        CoV = NaN;
    end

    % --- Eval count ---
    neval = eval_counter('get', 0);

    % --- Result struct ---
    result = create_result_struct();

    result.method    = 'MCS';
    result.Pf        = Pf;
    result.beta      = beta;
    result.neval     = neval;
    result.runtime   = toc;
    result.converged = true;
    result.COV = CoV;

    % --- Store diagnostics ---
    result.history.CoV = CoV;
    result.history.failures = sum(indicator);
    result.history.N = N;
    result.u_star = NaN;
    result.kappa  = NaN;

end