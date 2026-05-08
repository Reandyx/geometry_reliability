function b = calibrate_b(a, Pf_target, model_type)

    % ============================
    % DEFAULT MODEL TYPE
    % ============================
    if nargin < 3
        model_type = 'global'; % default
    end

    % ============================
    % INITIAL BRACKET
    % ============================
    b_low = 0;
    b_high = 5;

    N_mc_bracket = 1e5;

    % --- evaluate bounds ---
    problem_low = build_problem(a, b_low, model_type);
    Pf_low = run_mcs(problem_low, N_mc_bracket).Pf;

    problem_high = build_problem(a, b_high, model_type);
    Pf_high = run_mcs(problem_high, N_mc_bracket).Pf;

    % ============================
    % EXPAND BRACKET
    % ============================
    while ~(Pf_low > Pf_target && Pf_high < Pf_target)

        b_high = b_high * 2;

        problem_high = build_problem(a, b_high, model_type);
        Pf_high = run_mcs(problem_high, N_mc_bracket).Pf;

        % safety break
        if b_high > 1e3
            warning('Calibration failed: b_high too large');
            break;
        end
    end

    % ============================
    % BISECTION
    % ============================
    N_mc = 2e5;

    for k = 1:20

        b = 0.5 * (b_low + b_high);

        problem = build_problem(a, b, model_type);
        Pf = run_mcs(problem, N_mc).Pf;

        % DEBUG PRINT (keep this)
        % fprintf('model=%s | a=%.2f | iter=%d | b=%.4f | Pf=%.4f\n', ...
        %         model_type, a, k, b, Pf);

        if Pf > Pf_target
            b_low = b;
        else
            b_high = b;
        end

    end

end