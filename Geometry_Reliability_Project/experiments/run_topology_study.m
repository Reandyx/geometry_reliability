function results = run_topology_study()

    fprintf('\n============================\n');
    fprintf('   TOPOLOGY STUDY (M7)\n');
    fprintf('============================\n');

    target_pf = 0.15;
    N_mc = 1e6;

    % =========================
    % Calibration
    % =========================
    b = calibrate_b_disconnected(target_pf);

    fprintf('\nCalibrated b = %.4f\n', b);

    problem = get_problem_disconnected(b);

    % =========================
    % Monte Carlo
    % =========================
    reset_eval_count();
    res_mc = run_mcs(problem, N_mc);

    Pf_mc = res_mc.Pf;
    sigma_mc = sqrt(Pf_mc * (1 - Pf_mc) / N_mc);

    fprintf('\nMCS  Pf = %.6e ± %.2e\n', Pf_mc, sigma_mc);
    
    % =========================
    % FORM
    % =========================
    reset_eval_count();
    res_form = run_form(problem);

    Pf_form = res_form.Pf;

    fprintf('FORM Pf = %.6e\n', Pf_form);

    % =========================
    % SORM
    % =========================
    res_sorm = run_sorm(problem);

    Pf_sorm = res_sorm.Pf;

    fprintf('SORM Pf = %.6e\n', Pf_sorm);
    
    % =========================
    % Region + capture analysis
    % =========================
    N_split = 2e5;
    X = problem.sample(N_split);
    
    g_vals = evaluate_limit_state(problem, X);
    
    left_idx  = X(:,1) < 0;
    right_idx = X(:,1) >= 0;
    
    Pf_left  = mean((g_vals <= 0) & left_idx);
    Pf_right = mean((g_vals <= 0) & right_idx);
    
    dominant_region = max(Pf_left, Pf_right);
    capture_ratio = res_form.Pf / dominant_region;
    
    fprintf('Pf_left  = %.6e\n', Pf_left);
    fprintf('Pf_right = %.6e\n', Pf_right);
    fprintf('Sum      = %.6e\n', Pf_left + Pf_right);
    
    fprintf('Dominant region Pf ≈ %.6e\n', dominant_region);
    fprintf('FORM captures %.2f%% of dominant region\n', 100*capture_ratio);

    % =========================
    % Errors
    % =========================
    error_FORM = abs(Pf_form - Pf_mc) / Pf_mc;
    error_SORM = abs(Pf_sorm - Pf_mc) / Pf_mc;

    fprintf('\nError FORM = %.4f\n', error_FORM);
    fprintf('Error SORM = %.4f\n', error_SORM);
    
    missed_FORM = 1 - res_form.Pf / res_mc.Pf;
    missed_SORM = 1 - res_sorm.Pf / res_mc.Pf;
    
    fprintf('Missed FORM = %.4f\n', missed_FORM);
    fprintf('Missed SORM = %.4f\n', missed_SORM);

    % =========================
    % Topology insight
    % =========================
    U_star = res_form.U_star;

    if U_star(1) > 0
        region = 'RIGHT (+c)';
    else
        region = 'LEFT (-c)';
    end

    fprintf('\nFORM selected region: %s\n', region);
    fprintf('Second region ignored → topology failure\n');

    % =========================
    % Store
    % =========================
    results = struct();

    results.Pf_MC = Pf_mc;
    results.Pf_FORM = Pf_form;
    results.Pf_SORM = Pf_sorm;
    results.Pf_left  = Pf_left;
    results.Pf_right = Pf_right;

    results.error_FORM = error_FORM;
    results.error_SORM = error_SORM;
    results.missed_FORM = missed_FORM;
    results.missed_SORM = missed_SORM;

    results.U_star = U_star;
    results.region_selected = region;

    results.problem = problem;

    fprintf('\n============================\n');
    fprintf('   M7 COMPLETE\n');
    fprintf('============================\n');

end