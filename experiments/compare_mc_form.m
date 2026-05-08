function results_table = compare_mc_form(problem, N)
% COMPARE_MC_FORM
% Runs Monte Carlo and FORM and compares results

    % --- Run MCS ---
    res_mc = run_mcs(problem, N);

    % --- Run FORM ---
    res_form = run_form(problem);

    % --- Compute relative error ---
    if res_mc.Pf == 0
        error_rel = NaN;
    else
        error_rel = abs(res_form.Pf - res_mc.Pf) / res_mc.Pf;
    end

    % --- Store results ---
    results_table = struct();

    results_table.problem    = problem.name;
    results_table.Pf_MC      = res_mc.Pf;
    results_table.Pf_FORM    = res_form.Pf;
    results_table.error      = error_rel;
    results_table.beta_FORM  = res_form.beta;

    % --- Print results ---
    fprintf('\n============================\n');
    fprintf('   MCS vs FORM COMPARISON\n');
    fprintf('============================\n');

    fprintf('Problem: %s\n', problem.name);
    fprintf('Pf_MC   = %.5e\n', res_mc.Pf);
    fprintf('Pf_FORM = %.5e\n', res_form.Pf);
    fprintf('Error   = %.3e\n', error_rel);
    fprintf('beta_FORM = %.4f\n', res_form.beta);

end