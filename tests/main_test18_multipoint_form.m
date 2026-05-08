function main_test18_multipoint_form()

results = run_multipoint_study();

plot_multipoint_results(results);

fprintf('\n');
fprintf('=====================================\n');
fprintf('MULTI-POINT FORM SUMMARY\n');
fprintf('=====================================\n');

for i = 1:length(results)

    fprintf('\n');

    fprintf('c = %.2f\n', results(i).c);

    fprintf('Single capture = %.4f\n', ...
        results(i).capture_single);

    fprintf('Multi capture  = %.4f\n', ...
        results(i).capture_multi);

    fprintf('MCS Pf         = %.4e\n', ...
        results(i).MCS.Pf);

    fprintf('FORM Pf        = %.4e\n', ...
        results(i).FORM.Pf);

    fprintf('MULTI Pf       = %.4e\n', ...
        results(i).MULTI.Pf);

    fprintf('Bound lower    = %.4e\n', ...
        results(i).MULTI.bounds.lower);

    fprintf('Bound upper    = %.4e\n', ...
        results(i).MULTI.bounds.upper);

    fprintf('Branches       = %d\n', ...
        results(i).MULTI.n_branches);

end

fprintf('\n');
fprintf('A3 multipoint FORM study completed.\n');

end