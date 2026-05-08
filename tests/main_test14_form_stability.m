function main_test14_form_stability()

    fprintf('\n========================================\n');
    fprintf(' TEST 14 — FORM STABILITY\n');
    fprintf('========================================\n');

    n_trials = 20;
    converged_count = 0;
    residuals = zeros(n_trials,1);
    betas = zeros(n_trials,1);

    % Base problem
    base_problem = get_problem_synthetic();

    for i = 1:n_trials

        % Random initial guess (stress test)
        problem = base_problem;
        problem.U0 = randn(problem.dimension,1) * 2;

        res = run_form(problem);

        residuals(i) = res.residual;
        betas(i) = res.beta;

        if res.converged && res.residual < 1e-6
            converged_count = converged_count + 1;
        end

        fprintf('Trial %d | converged=%d | residual=%.2e | beta=%.4f\n', ...
            i, res.converged, res.residual, res.beta);

    end

    convergence_rate = converged_count / n_trials;

    fprintf('\n--- SUMMARY ---\n');
    fprintf('Convergence rate = %.2f\n', convergence_rate);
    fprintf('Max residual     = %.2e\n', max(residuals));
    fprintf('Std(beta)        = %.4e\n', std(betas));

    % ==========================
    % PASS / FAIL CONDITIONS
    % ==========================
    assert(convergence_rate >= 0.95, 'FAIL: convergence rate too low');
    assert(max(residuals) < 1e-5, 'FAIL: residual too large');
    assert(std(betas) < 1e-2, 'FAIL: solution depends on initial guess');

    fprintf('\n✅ FORM STABILITY TEST PASSED\n');

end