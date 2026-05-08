function results = run_distribution_sensitivity()

rng(1);

% =====================================
% Base problem parameters
% =====================================
a = 2.0;
b = 4.0;

N_mcs = 1e5;

% =====================================
% Distribution sets
% =====================================
dist_sets = {
    {'normal','normal'}
    {'lognormal','lognormal'}
    {'gumbel','gumbel'}
};

dist_names = {
    'normal'
    'lognormal'
    'gumbel'
};

results = struct();

for i = 1:length(dist_sets)

    % =====================================
    % Problem
    % =====================================
    problem = get_problem_local_curvature(a, b);

    % IMPORTANT:
    % lognormal requires positive mean
    if strcmp(dist_names{i}, 'lognormal')

        problem.mu = [5 5];
        problem.sigma = [1 1];

    else

        problem.mu = [0 0];
        problem.sigma = [1 1];

    end

    problem.dist = dist_sets{i};

    % =====================================
    % Methods
    % =====================================
    res_form = run_form(problem);

    res_sorm = run_sorm(problem);

    res_mcs = run_mcs(problem, N_mcs);

    % =====================================
    % Geometry metrics
    % =====================================
    metrics = compute_curvature_metrics( ...
        res_form.beta, ...
        res_sorm.kappa);

    % =====================================
    % Store
    % =====================================
    results(i).distribution = dist_names{i};

    results(i).FORM = res_form;
    results(i).SORM = res_sorm;
    results(i).MCS  = res_mcs;

    results(i).gamma = metrics.gamma;

    results(i).beta = res_form.beta;

    results(i).Pf = res_form.Pf;

end

% =====================================
% Beta shift relative to normal
% =====================================
beta_normal = results(1).beta;

for i = 1:length(results)

    results(i).delta_beta = ...
        results(i).beta - beta_normal;

    results(i).Pf_ratio = ...
        results(i).Pf / results(1).Pf;

end

% =====================================
% Save
% =====================================
if ~exist('results/distribution/raw','dir')
    mkdir('results/distribution/raw');
end

save('results/distribution/raw/data.mat','results');

end