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

    % --- CoV ---
    if Pf == 0
        CoV = Inf;
    else
        var_est = var(indicator) / N;
        CoV = sqrt(var_est) / Pf;
    end

    % --- Reliability index ---
    if Pf == 0
        beta = Inf;
    else
        beta = compute_beta_from_pf(Pf);
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

    % --- Store diagnostics ---
    result.history.CoV = CoV;
    result.history.failures = sum(indicator);
    result.history.N = N;

end