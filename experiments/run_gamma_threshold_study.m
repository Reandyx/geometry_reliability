function results = run_gamma_threshold_study()

rng(1);

% =====================================
% Curvature sweep
% =====================================
a_values = [0.25 0.5 1 2 3 5 7 10 15 20];

Pf_target = 0.20;

N_mcs = 1e5;

results = struct();

for i = 1:length(a_values)

    a = a_values(i);

    % =====================================
    % Benchmark
    % =====================================
    b = calibrate_b(a, Pf_target, 'local');

    problem = get_problem_local_curvature(a, b);

    % =====================================
    % Methods
    % =====================================
    res_form = run_form(problem);

    res_sorm = run_sorm(problem);

    res_mcs = run_mcs(problem, N_mcs);

    % =====================================
    % Geometry index
    % =====================================
    gamma = compute_geometry_index( ...
        res_form.beta, ...
        res_sorm.kappa);

    % =====================================
    % Errors
    % =====================================
    form_error = ...
        abs(res_form.Pf - res_mcs.Pf) ...
        / res_mcs.Pf;

    sorm_error = ...
        abs(res_sorm.Pf - res_mcs.Pf) ...
        / res_mcs.Pf;

    % =====================================
    % Improvement ratio
    % =====================================
    if form_error > 0

        sorm_gain = ...
            (form_error - sorm_error) ...
            / form_error;

    else

        sorm_gain = 0;

    end

    % =====================================
    % Classification
    % =====================================
    if gamma < 0.1
        regime = 'Low-curvature';
    elseif gamma < 1
        regime = 'Moderate-curvature';
    else
        regime = 'High-curvature';
    end

    % =====================================
    % Store
    % =====================================
    results(i).a = a;

    results(i).gamma = gamma;

    results(i).FORM = res_form;

    results(i).SORM = res_sorm;

    results(i).MCS = res_mcs;

    results(i).form_error = form_error;

    results(i).sorm_error = sorm_error;

    results(i).sorm_gain = sorm_gain;

    results(i).regime = regime;

end

% =====================================
% Save
% =====================================
if ~exist('results/gamma_thresholds/raw','dir')

    mkdir('results/gamma_thresholds/raw');

end

save('results/gamma_thresholds/raw/data.mat', ...
    'results');

end