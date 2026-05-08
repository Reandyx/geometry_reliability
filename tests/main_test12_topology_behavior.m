function main_test12_topology_behavior()

    fprintf('\n========================================\n');
    fprintf(' TEST 12B — TOPOLOGY BEHAVIOR\n');
    fprintf('========================================\n');

    c = 2.0;
    target_pf = 0.15;

    b = calibrate_b_disconnected_c(c, target_pf);
    problem = get_problem_disconnected_c(b, c);

    % ==========================
    % MCS (truth)
    % ==========================
    res_mc = run_mcs(problem, 5e5);

    % ==========================
    % FORM
    % ==========================
    res_form = run_form(problem);

    ratio = res_form.Pf / res_mc.Pf;

    fprintf('Pf_MC   = %.6e\n', res_mc.Pf);
    fprintf('Pf_FORM = %.6e\n', res_form.Pf);
    fprintf('Capture = %.2f%%\n', 100*ratio);

    % Expected: FORM underestimates
    assert(ratio < 1.0, 'FAIL: FORM should underestimate');

    % But not collapse completely
    assert(ratio > 0.3, 'FAIL: FORM too inaccurate');

    fprintf('\n✅ TEST PASSED\n');

end