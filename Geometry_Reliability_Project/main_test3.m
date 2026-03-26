clc;
clear;

addpath(genpath(pwd));
rng(1);

fprintf('============================\n');
fprintf('        M3 Test\n');
fprintf('============================\n\n');

%% --- Cantilever Test ---
problem = get_problem_cantilever();

N = 1e5;

res1 = compare_mc_form(problem, N);

%% --- Synthetic Curvature Study ---
results = [];

a_values = [0.5, 1.0, 2.0, 5.0];

for i = 1:length(a_values)

    problem = get_problem_synthetic(a_values(i));

    res = compare_mc_form(problem, 1e5);

    res.a = a_values(i);
    res.b = problem.b;

    results = [results; res];

end