function problem = get_problem_cantilever()

    % Name
    problem.name = 'Cantilever Beam';

    % Dimension: X = [E, F]
    problem.dimension = 2;

    % Mean values
    problem.mu = [210e9, 1000];   % [E, F]

    % Standard deviations
    problem.sigma = [21e9, 200];

    % Distribution type (assume normal for now)
    problem.dist = {'normal', 'normal'};

    % Deterministic parameters
    L = 2.0;
    b = 0.05;
    h = 0.10;
    delta_allow = 0.004;

    I = b * h^3 / 12;

    % Limit-state function g(X)
    problem.gfun = @(X) cantilever_gfun(X, L, I, delta_allow);

    % Sampling handle
    problem.sample = @(N) sample_inputs(problem, N);

    problem.is_normal = true;
    problem.is_independent = true;

end

function g = cantilever_gfun(X, L, I, delta_allow)

    E = X(:,1);
    F = X(:,2);

    g = delta_allow - (F .* L^3) ./ (3 .* E .* I);

end