function b_opt = calibrate_b_disconnected(target_pf)

    if nargin < 1
        target_pf = 0.15;
    end

    b_low = 0.5;
    b_high = 6.0;

    tol = 1e-3;
    max_iter = 20;

    for i = 1:max_iter

        b_mid = 0.5*(b_low + b_high);

        problem = get_problem_disconnected(b_mid);

        res = run_mcs(problem, 2e5); % lighter MC for calibration

        Pf = res.Pf;

        fprintf('[CAL] iter=%d | b=%.4f | Pf=%.4f\n', i, b_mid, Pf);

        if abs(Pf - target_pf) < tol
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