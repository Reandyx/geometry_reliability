function b_opt = calibrate_b_disconnected_c(c, target_pf)

    % ============================
    % INITIAL GUESS
    % ============================
    b_low = 1.0;
    b_high = 2.0;

    N_mc_bracket = 1e5;

    % ============================
    % FIND VALID BRACKET
    % ============================
    bracket_found = false;

    for iter = 1:25

        problem_low = get_problem_disconnected_c(b_low, c);
        Pf_low = run_mcs(problem_low, N_mc_bracket).Pf;

        problem_high = get_problem_disconnected_c(b_high, c);
        Pf_high = run_mcs(problem_high, N_mc_bracket).Pf;

        % CORRECT CONDITION
        if Pf_low < target_pf && Pf_high > target_pf
            bracket_found = true;
            break;
        end

        % expand search (both directions)
        b_low = b_low / 2;
        b_high = b_high * 2;

        if b_high > 1e6
            break;
        end

    end

    % ============================
    % HANDLE FAILURE
    % ============================
    if ~bracket_found
        warning('Calibration impossible for c=%.2f (Pf never crosses target)', c);
        b_opt = NaN;
        return;
    end

    % ============================
    % BISECTION
    % ============================
    for i = 1:25

        b_mid = 0.5 * (b_low + b_high);

        problem = get_problem_disconnected_c(b_mid, c);
        Pf = run_mcs(problem, 2e5).Pf;

        if abs(Pf - target_pf) < 1e-3
            b_opt = b_mid;
            return;
        end

        if Pf > target_pf
            b_high = b_mid;
        else
            b_low = b_mid;
        end

    end

    b_opt = b_mid;

end