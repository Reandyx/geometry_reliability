function problem = get_problem_truss()

problem.name = 'Two-bar truss';
problem.model_type = 'engineering';

problem.dimension = 2;

% Random variables: cross-sectional areas
problem.mu = [0.01, 0.01];
problem.sigma = [0.002, 0.002];
problem.dist = {'normal','normal'};

% Deterministic load
P = 2;

problem.gfun = @(X) truss_gfun(X, P);
problem.sample = @(N) sample_inputs(problem, N);

problem.is_normal = true;
problem.is_independent = true;

end


function g = truss_gfun(X, P)

A1 = X(:,1);
A2 = X(:,2);

stress = P ./ (A1 + A2);

sigma_allow = 250;

g = sigma_allow - stress;

end