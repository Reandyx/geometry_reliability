function main_test13_curvature_consistency()

    fprintf('\n========================================\n');
    fprintf(' TEST 13 — CURVATURE CONSISTENCY\n');
    fprintf('========================================\n');

    problem = get_problem_synthetic();

    % ==========================
    % FORM
    % ==========================
    res_form = run_form(problem);

    beta = res_form.beta;
    U_star = res_form.U_star;
    grad_g = res_form.grad_g;

    % ==========================
    % GEOMETRY PIPELINE
    % ==========================
    H = compute_hessian_u(problem, U_star);

    kappa = compute_principal_curvatures(H, grad_g);

    metrics = compute_curvature_metrics(beta, kappa);

    % ==========================
    % SORM (internal)
    % ==========================
    res_sorm = run_sorm(problem);

    Pf_form = res_form.Pf;
    Pf_sorm = res_sorm.Pf;

    % ==========================
    % MANUAL BREITUNG
    % ==========================
    prod_term = 1;

    for i = 1:length(kappa)
        val = 1 + beta * kappa(i);
        assert(val > 0, 'SORM breakdown in manual check');
        prod_term = prod_term * val^(-0.5);
    end

    Pf_manual = Pf_form * prod_term;

    % ==========================
    % CHECK CONSISTENCY
    % ==========================
    rel_error = abs(Pf_manual - Pf_sorm) / Pf_sorm;

    fprintf('Pf_FORM   = %.6e\n', Pf_form);
    fprintf('Pf_SORM   = %.6e\n', Pf_sorm);
    fprintf('Pf_manual = %.6e\n', Pf_manual);
    fprintf('Rel error = %.3e\n', rel_error);

    fprintf('\n--- CURVATURE METRICS ---\n');
    fprintf('kappa_max = %.6f\n', metrics.kappa_max);
    fprintf('gamma     = %.6f\n', metrics.gamma);

    % ==========================
    % PASS / FAIL
    % ==========================
    assert(rel_error < 1e-6, 'FAIL: SORM mismatch with manual Breitung');
    assert(metrics.gamma >= 0, 'FAIL: gamma must be non-negative');

    fprintf('\n✅ CURVATURE CONSISTENCY TEST PASSED\n');

end