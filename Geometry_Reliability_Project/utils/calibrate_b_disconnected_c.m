function b_opt = calibrate_b_disconnected_c(c, target_pf)

    b_low = 0.5;
    b_high = 12.0;

    for i = 1:20

        b_mid = 0.5*(b_low + b_high);

        problem = get_problem_disconnected_c(b_mid, c);

        res = run_mcs(problem, 2e5);
        Pf = res.Pf;

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