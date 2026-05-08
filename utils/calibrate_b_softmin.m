function b = calibrate_b_softmin(c, Pf_target, k)

    b_low = 0;
    b_high = 10;

    for iter = 1:25

        b_mid = 0.5*(b_low + b_high);

        beta1 = b_mid;
        beta2 = b_mid + 1.0;

        %SAME MODEL AS EXPERIMENT (SOFTMIN)
        problem = get_problem_softmin(beta1, beta2, c, k);

        res_mcs = run_mcs(problem, 2e5);

        if res_mcs.Pf > Pf_target
            b_low = b_mid;
        else
            b_high = b_mid;
        end

    end

    b = b_mid;

end