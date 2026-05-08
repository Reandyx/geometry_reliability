function bounds = compute_ditlevsen_bounds(branch_results)

n = length(branch_results);

Pf = zeros(n,1);

for i = 1:n
    Pf(i) = branch_results(i).Pf;
end

% =====================================
% Simple bounds
% =====================================

Pf_upper = sum(Pf);

Pf_lower = max(Pf);

bounds.lower = Pf_lower;
bounds.upper = Pf_upper;

end