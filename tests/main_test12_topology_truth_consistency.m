function main_test12_topology_truth_consistency()
    
    fprintf('\n========================================\n');
    fprintf(' TEST 12A — TOPOLOGY TRUTH CONSISTENCY\n');
    fprintf('========================================\n');

    c = 2.0;
    target_pf = 0.15;

    b = calibrate_b_disconnected_c(c, target_pf);
    problem = get_problem_disconnected_c(b, c);

    % ==========================
    % Random points
    % ==========================
    U = randn(1000,2);

    % Direct evaluation
    g_direct = problem.gfun(U);

    % Manual reconstruction
    a = problem.a;
    b_val = problem.b;
    c_val = problem.c;
    k = problem.k;

    g1 = (U(:,2).^2 + a*(U(:,1) - c_val).^2) - b_val;
    g2 = (U(:,2).^2 + a*(U(:,1) + c_val).^2) - b_val + 0.05;

    m = min([g1 g2],[],2);

    g_manual = m - (1/k)*log(exp(-k*(g1-m)) + exp(-k*(g2-m)));

    err = max(abs(g_direct - g_manual));

    fprintf('Max g mismatch = %.3e\n', err);

    assert(err < 1e-10, 'FAIL: inconsistent g definition');

    fprintf('\n✅ TEST PASSED\n');

end