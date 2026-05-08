function main_test18_distribution_sensitivity()

results = run_distribution_sensitivity();

plot_distribution_sensitivity(results);

fprintf('\n');
fprintf('=====================================\n');
fprintf('DISTRIBUTION SENSITIVITY SUMMARY\n');
fprintf('=====================================\n');

for i = 1:length(results)

    fprintf('\n');

    fprintf('Distribution : %s\n', ...
        results(i).distribution);

    fprintf('beta         : %.4f\n', ...
        results(i).beta);

    fprintf('delta beta   : %.4f\n', ...
        results(i).delta_beta);

    fprintf('Pf           : %.4e\n', ...
        results(i).Pf);

    fprintf('Pf ratio     : %.4f\n', ...
        results(i).Pf_ratio);

    fprintf('gamma        : %.4f\n', ...
        results(i).gamma);

    fprintf('iterations   : %d\n', ...
        results(i).FORM.n_iterations);

    fprintf('cond(H)      : %.4e\n', ...
        results(i).FORM.cond_H);

    fprintf('residual     : %.4e\n', ...
        results(i).FORM.residual);

    fprintf('converged    : %d\n', ...
        results(i).FORM.converged);

end

fprintf('\n');
fprintf('A2 distribution sensitivity study completed.\n');

end