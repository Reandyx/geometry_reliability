function problem = get_problem_local_curvature(a, b)

    if nargin < 1
        a = 1.0;
    end
    if nargin < 2
        b = 4.0;   % threshold (controls Pf)
    end

    problem.name = sprintf('Synthetic Local Curvature (a=%.2f)', a);
    problem.model_type = 'local';

    problem.dimension = 2;

    problem.mu = [0, 0];
    problem.sigma = [1, 1];
    problem.dist = {'normal', 'normal'};

    problem.a = a;
    problem.b = b;   

    problem.gfun = @(X) local_gfun(X, problem);
    problem.sample = @(N) sample_inputs(problem, N);

    problem.is_normal = true;
    problem.is_independent = true;

end


function g = local_gfun(X, problem)

    x1 = X(:,1);
    x2 = X(:,2);

    a = problem.a;
    b = problem.b;

    g = b - (x1.^2 + a * exp(-x1.^2) .* x2.^2 + 0.1*x1);

end