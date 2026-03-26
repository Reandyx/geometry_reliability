function problem = get_problem_synthetic(a, b)

    if nargin < 1
        a = 1.0;
    end

    if nargin < 2
        b = 4.0; % threshold (controls Pf)
    end

    problem.name = sprintf('Synthetic Global Curvature (a=%.2f)', a);
    problem.model_type = 'global';
    problem.dimension = 2;

    problem.mu = [0, 0];
    problem.sigma = [1, 1];

    problem.dist = {'normal', 'normal'};

    problem.a = a;  
    problem.b = b;  

    problem.gfun = @(X) synthetic_gfun(X, problem);

    problem.sample = @(N) sample_inputs(problem, N);

    % FLAGS
    problem.is_normal = true;
    problem.is_independent = true;

end

function g = synthetic_gfun(X, problem)

    x1 = X(:,1);
    x2 = X(:,2);

    a = problem.a;
    b = problem.b;

    g = b - ( (x1 - 1).^2 + a * (x2 - 0.5).^2 );

end