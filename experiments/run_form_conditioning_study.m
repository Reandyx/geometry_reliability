function results = run_form_conditioning_study()

rng(1);

% =====================================
% Curvature sweep
% =====================================
a_values = [0.25 0.5 1 2 3 5 7 10 15 20];

Pf_target = 0.2;
N_mcs = 1e5;

results = struct();

for i = 1:length(a_values)

    a = a_values(i);

    % --- calibration ---
    b = calibrate_b(a, Pf_target, 'local');

    % --- problem ---
    problem = build_problem(a, b, 'local');

    % --- FORM ---
    res_form = run_form(problem);

    % --- SORM ---
    res_sorm = run_sorm(problem);

    % --- MCS ---
    res_mcs = run_mcs(problem, N_mcs);

    % --- geometry metrics ---
    metrics = compute_curvature_metrics( ...
        res_form.beta, ...
        res_sorm.kappa);

    % =====================================
    % Error metrics
    % =====================================
    form_error = abs(res_form.Pf - res_mcs.Pf) / res_mcs.Pf;

    if ~isnan(res_sorm.Pf)
        sorm_error = abs(res_sorm.Pf - res_mcs.Pf) / res_mcs.Pf;
    else
        sorm_error = NaN;
    end

    % =====================================
    % Store
    % =====================================
    results(i).a = a;

    results(i).gamma = metrics.gamma;

    results(i).FORM = res_form;
    results(i).SORM = res_sorm;
    results(i).MCS  = res_mcs;

    results(i).form_error = form_error;
    results(i).sorm_error = sorm_error;

    results(i).cond_H = res_form.cond_H;
    results(i).iterations = res_form.n_iterations;
    results(i).residual = res_form.residual;

end

% =====================================
% Save
% =====================================
if ~exist('results/conditioning/raw', 'dir')
    mkdir('results/conditioning/raw');
end

save('results/conditioning/raw/data.mat', 'results');

fprintf('\n');
fprintf('=====================================\n');
fprintf('FORM CONDITIONING SUMMARY\n');
fprintf('=====================================\n');

for i = 1:length(results)

    fprintf(['a=%6.2f | gamma=%6.3f | cond(H)=%10.3e | ' ...
             'iter=%3d | FORM err=%6.3f | SORM err=%6.3f\n'], ...
             results(i).a, ...
             results(i).gamma, ...
             results(i).cond_H, ...
             results(i).iterations, ...
             results(i).form_error, ...
             results(i).sorm_error);

end

end