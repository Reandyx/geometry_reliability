function results = run_multipoint_study()

rng(1);

% =====================================
% Topology parameter sweep
% =====================================
c_values = [0.5 1.0 2.0 3.0];

target_pf = 0.15;

N_mcs = 1e5;

results = struct();

for i = 1:length(c_values)

    c = c_values(i);

    % =====================================
    % Calibrate benchmark
    % =====================================
    b = calibrate_b_disconnected_c(c, target_pf);

    problem = get_problem_disconnected_c(b, c);

    % =====================================
    % Single FORM
    % =====================================
    res_form = run_form(problem);

    % =====================================
    % Multi-point FORM
    % =====================================
    res_multi = run_multipoint_form(problem);

    % =====================================
    % MCS reference
    % =====================================
    res_mcs = run_mcs(problem, N_mcs);

    % =====================================
    % Capture ratios
    % =====================================
    capture_single = ...
        res_form.Pf / res_mcs.Pf;

    capture_multi = ...
        res_multi.Pf / res_mcs.Pf;

    % =====================================
    % Store
    % =====================================
    results(i).c = c;

    results(i).FORM = res_form;

    results(i).MULTI = res_multi;

    results(i).MCS = res_mcs;

    results(i).capture_single = capture_single;

    results(i).capture_multi = capture_multi;

end

% =====================================
% Save
% =====================================
if ~exist('results/multipoint/raw','dir')
    mkdir('results/multipoint/raw');
end

save('results/multipoint/raw/data.mat','results');

end