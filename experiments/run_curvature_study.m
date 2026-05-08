function results = run_curvature_study()

rng(1);

% ========================
% CONFIG
% ========================
a_values = [0.25 0.5 1 2 5 10 20 50]
Pf_target = 0.2;
N_mcs = 1e5;

results = struct();

for i = 1:length(a_values)

a = a_values(i);

% --- calibration ---
b = calibrate_b(a, Pf_target, 'global');

% --- problem ---
problem = build_problem(a, b, 'global');

% --- methods ---
res_form = run_form(problem);
res_sorm = run_sorm(problem);
res_mcs  = run_mcs(problem, N_mcs);

% --- curvature ---
metrics = compute_curvature_metrics(res_form.beta, res_sorm.kappa);

% --- store ---
results(i).a = a;
results(i).gamma = metrics.gamma;

results(i).FORM = res_form;
results(i).SORM = res_sorm;
results(i).MCS  = res_mcs;

end

% --- save ---
if ~exist('results/curvature/raw', 'dir')
mkdir('results/curvature/raw');
end

save('results/curvature/raw/data.mat', 'results');

end
