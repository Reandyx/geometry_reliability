function results = run_benchmark_cases()

rng(1);

cases = {
    'cantilever', @get_problem_cantilever;
    'truss',      @get_problem_truss;
    'nonlinear',  @get_problem_nonlinear;
    'borderline', @() get_problem_synthetic(2.0, 4, 0.0);
    'high_curv',  @() get_problem_synthetic(5.0, 4, 0.0);
};

for i = 1:size(cases,1)

    name = cases{i,1};
    problem_fun = cases{i,2};

    problem = problem_fun();

    % --- methods ---
    res_form = run_form(problem);
    res_sorm = run_sorm(problem);
    res_mcs  = run_mcs(problem, 1e6);

    % --- curvature ---
    metrics = compute_curvature_metrics(res_form.beta, res_sorm.kappa);
    gamma = metrics.gamma;

    % --- topology (default: connected) ---
    topology_flag = 0;

    % --- predicted method ---
    class = classify_geometry(gamma, topology_flag, res_mcs.Pf);
    
    predicted = select_method_from_geometry( ...
        class, ...
        gamma, ...
        res_form.beta, ...
        topology_flag ...
    );

    % --- errors ---
    err_form = abs(res_form.Pf - res_mcs.Pf) / res_mcs.Pf;
    err_sorm = abs(res_sorm.Pf - res_mcs.Pf) / res_mcs.Pf;

    % --- best method ---
    if res_form.beta >= 5
        best = 'IS';
    elseif abs(err_form - err_sorm) < 0.01
        best = 'FORM';   % prefer simpler method
    elseif err_form < err_sorm
        best = 'FORM';
    else
        best = 'SORM';
    end

    % --- store ---
    results(i).name = name;
    results(i).gamma = gamma;
    results(i).beta = res_form.beta;
    results(i).FORM = res_form;
    results(i).SORM = res_sorm;
    results(i).MCS  = res_mcs;
    results(i).err_form = err_form;
    results(i).err_sorm = err_sorm;
    results(i).predicted = predicted;
    results(i).best = best;

end

if ~exist('results/benchmark/raw','dir')
    mkdir('results/benchmark/raw');
end

save('results/benchmark/raw/data.mat','results');

end