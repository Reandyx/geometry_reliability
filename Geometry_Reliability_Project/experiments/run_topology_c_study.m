function results = run_topology_c_study()

    fprintf('\n============================\n');
    fprintf('   TOPOLOGY TRANSITION STUDY (c-variation)\n');
    fprintf('============================\n');

    % Sweep values
    c_values = [0.5, 1.0, 1.5, 2.0, 3.0];
    target_pf = 0.15;
    N_mc = 1e6;

    n_cases = length(c_values);

    results(n_cases) = struct( ...
        'c', [], ...
        'Pf_MC', [], ...
        'Pf_FORM', [], ...
        'Pf_SORM', [], ...
        'error_FORM', [], ...
        'error_SORM', [], ...
        'beta_FORM', [], ...
        'U_star', [], ...
        'region_selected', [], ...
        'missed_FORM', [], ...
        'missed_SORM', [], ...
        'Pf_left', [], ...
        'Pf_right', [], ...
        'dominant_region', [], ...
        'capture_ratio', []);

    for i = 1:n_cases

        c = c_values(i);

        fprintf('\n--- CASE c = %.2f ---\n', c);

        % =========================
        % Calibrate b
        % =========================
        b = calibrate_b_disconnected_c(c, target_pf);
        fprintf('Calibrated b = %.4f\n', b);

        problem = get_problem_disconnected_c(b, c);

        % =========================
        % Monte Carlo
        % =========================
        res_mc = run_mcs(problem, N_mc);
        Pf_mc = res_mc.Pf;

        fprintf('MCS  Pf = %.6e\n', Pf_mc);

        % =========================
        % Region contribution split (CORRECT)
        % =========================
        N_split = 2e5;
        X = problem.sample(N_split);
        g_vals = evaluate_limit_state(problem, X);

        left_idx  = X(:,1) < 0;
        right_idx = X(:,1) >= 0;

        Pf_left  = mean((g_vals <= 0) & left_idx);
        Pf_right = mean((g_vals <= 0) & right_idx);

        % =========================
        % FORM
        % =========================
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
        % Errors
        % =========================
        error_FORM = abs(Pf_form - Pf_mc) / Pf_mc;
        error_SORM = abs(Pf_sorm - Pf_mc) / Pf_mc;

        missed_FORM = 1 - Pf_form / Pf_mc;
        missed_SORM = 1 - Pf_sorm / Pf_mc;

        % =========================
        % Dominant region + capture
        % =========================
        dominant_region = max(Pf_left, Pf_right);
        capture_ratio = Pf_form / dominant_region;

        % =========================
        % Region detection
        % =========================
        U_star = res_form.U_star;

        if U_star(1) > 0
            region = 'RIGHT';
        else
            region = 'LEFT';
        end

        % =========================
        % Print diagnostics
        % =========================
        fprintf('Pf_left  = %.6e\n', Pf_left);
        fprintf('Pf_right = %.6e\n', Pf_right);
        fprintf('Sum      = %.6e\n', Pf_left + Pf_right);

        fprintf('Error FORM = %.4f\n', error_FORM);
        fprintf('Error SORM = %.4f\n', error_SORM);

        fprintf('Missed FORM = %.4f\n', missed_FORM);
        fprintf('Missed SORM = %.4f\n', missed_SORM);

        fprintf('Dominant region ≈ %.6e\n', dominant_region);
        fprintf('FORM captures %.2f%% of dominant region\n', 100*capture_ratio);

        % =========================
        % Store
        % =========================
        results(i).c = c;

        results(i).Pf_MC = Pf_mc;
        results(i).Pf_FORM = Pf_form;
        results(i).Pf_SORM = Pf_sorm;

        results(i).error_FORM = error_FORM;
        results(i).error_SORM = error_SORM;

        results(i).beta_FORM = res_form.beta;

        results(i).U_star = U_star;
        results(i).region_selected = region;

        results(i).missed_FORM = missed_FORM;
        results(i).missed_SORM = missed_SORM;

        results(i).Pf_left  = Pf_left;
        results(i).Pf_right = Pf_right;

        results(i).dominant_region = dominant_region;
        results(i).capture_ratio = capture_ratio;
        
        results(i).problem = problem;

    end

    fprintf('\n=== TOPOLOGY STUDY COMPLETE ===\n');

end