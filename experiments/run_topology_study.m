function results = run_topology_study()

rng(1);

c_values = [0.5, 1.0, 2.0, 3.0];
target_pf = 0.15;
N_mcs = 1e6;

results = struct();

for i = 1:length(c_values)

c = c_values(i);

% --- calibration ---
b = calibrate_b_disconnected_c(c, target_pf);

% --- problem ---
problem = get_problem_disconnected_c(b, c);

% --- methods ---
res_form = run_form(problem);
res_sorm = run_sorm(problem);
res_mcs  = run_mcs(problem, N_mcs);

% --- capture ratio ---
capture_ratio = res_form.Pf / res_mcs.Pf;

% --- store ---
results(i).c = c;
results(i).FORM = res_form;
results(i).SORM = res_sorm;
results(i).MCS  = res_mcs;
results(i).capture_ratio = capture_ratio;

end

if ~exist('results/topology/raw', 'dir')
mkdir('results/topology/raw');
end

save('results/topology/raw/data.mat', 'results');

end
