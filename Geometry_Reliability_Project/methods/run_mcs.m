function result = run_mcs(problem, N)
% RUN_MCS
% Crude Monte Carlo Simulation for reliability estimation

    tic;

    % Reset evaluation counter
    reset_eval_count();

    % 1. Generate samples (use problem interface)
    X = sample_inputs(problem, N);

    % 2. Evaluate limit-state (USE CORE PIPELINE)
    g_vals = evaluate_limit_state(problem, X);

    % 3. Failure probability
    Pf = mean(g_vals <= 0);

    % 4. Reliability index
    if Pf == 0
        beta = Inf;
    else
        beta = compute_beta_from_pf(Pf);
    end

    % 5. Number of evaluations
    neval = eval_counter('get', 0);

    % 6. Build result struct
    result = create_result_struct();

    result.method    = 'MCS';
    result.Pf        = Pf;
    result.beta      = beta;
    result.neval     = neval;
    result.runtime   = toc;
    result.converged = true;

    % Safety checks (interface enforcement)
    assert(isfield(result, 'Pf'));
    assert(isfield(result, 'beta'));
    assert(isfield(result, 'U_star'));
    
end