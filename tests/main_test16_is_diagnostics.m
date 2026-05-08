function main_test16_is_diagnostics()

    fprintf('\n========================================\n');
    fprintf(' TEST 16 — IS DIAGNOSTICS\n');
    fprintf('========================================\n');

    problem = get_problem_local_curvature(1.0, 10);

    res_form = run_form(problem);

    N = 1e5;

    % --- NEW INTERFACE ---
    res_is = run_is(problem, N, res_form.U_star);

    Pf_IS  = res_is.Pf;
    CoV_IS = res_is.history.CoV;
    ESS_ratio = res_is.history.ESS_ratio;

    fprintf('\nPf_IS = %.3e\n', Pf_IS);
    fprintf('CoV   = %.3f\n', CoV_IS);
    fprintf('ESS%%  = %.2f%%\n', 100*ESS_ratio);

    % =========================
    % SANITY CHECKS
    % =========================
    assert(~isnan(Pf_IS), 'FAIL: IS returned NaN');
    assert(ESS_ratio > 0.01, 'FAIL: severe degeneracy');

    fprintf('\n✅ IS TEST PASSED\n');

end