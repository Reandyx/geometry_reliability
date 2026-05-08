function result = run_multipoint_form(problem)

% =====================================
% Seed generation
% =====================================
c = problem.c;

seeds = {
    [ c;  0.5]
    [-c; 0.5]
};

branch_results = struct([]);

for i = 1:length(seeds)

    % =====================================
    % Set FORM initialization
    % =====================================
    problem_i = problem;

    problem_i.U0 = seeds{i};

    % =====================================
    % Run FORM
    % =====================================
    res = run_form(problem_i);

    if ~res.converged
        continue;
    end

    branch_results = [branch_results; res];

end

if isempty(branch_results)

    warning('Multipoint FORM: no converged branches');

    result = create_result_struct();

    result.method = 'MULTIPOINT_FORM';

    result.Pf = NaN;

    return;

end

% =====================================
% Duplicate filtering
% =====================================
unique_branches = struct([]);

tol_duplicate = 1e-2;

for i = 1:length(branch_results)

    is_duplicate = false;

    for j = 1:length(unique_branches)

        d = norm( ...
            branch_results(i).U_star - ...
            unique_branches(j).U_star);

        if d < tol_duplicate

            is_duplicate = true;
            break;

        end

    end

    if ~is_duplicate

        unique_branches = ...
            [unique_branches; branch_results(i)];

    end

end

branch_results = unique_branches;

% Ditlevsen-style bounds
bounds = compute_ditlevsen_bounds(branch_results);

% =====================================
% Aggregate probabilities
% =====================================
Pf_total = 0;

for i = 1:length(branch_results)

    Pf_total = Pf_total + branch_results(i).Pf;

end

% =====================================
% Build result
% =====================================
result = create_result_struct();

result.method = 'MULTIPOINT_FORM';
result.branch_results = branch_results;
result.Pf = Pf_total;
result.beta = compute_beta_from_pf(Pf_total);
result.n_branches = length(branch_results);
result.bounds = bounds;

end