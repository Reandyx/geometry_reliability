function main_test15_calibration_accuracy()

    fprintf('\n========================================\n');
    fprintf(' TEST 15 — CALIBRATION ACCURACY\n');
    fprintf('========================================\n');

    c_values = [0.5, 1.0, 2.0, 3.0];
    target_pf = 0.15;

    for i = 1:length(c_values)

        c = c_values(i);

        b = calibrate_b_disconnected_c(c, target_pf);
        problem = get_problem_disconnected_c(b, c);

        res = run_mcs(problem, 5e5);

        fprintf('c=%.2f | target=%.3f | actual=%.3f\n', ...
            c, target_pf, res.Pf);

    end

end