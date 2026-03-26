function problem = get_problem_disconnected(b)

    if nargin < 1
        b = 2.5; % initial guess (will be calibrated)
    end

    a = 1.5;
    c = 2.0;

    problem.name = 'Disconnected Parabolic Domain';
    problem.model_type = 'topology';
    problem.dimension = 2;

    problem.mu = [0, 0];
    problem.sigma = [1, 1];

    problem.dist = {'normal','normal'};

    problem.a = a;
    problem.b = b;
    problem.c = c;

    problem.gfun = @(X) disconnected_gfun(X, problem);

    problem.sample = @(N) sample_inputs(problem, N);

    problem.is_normal = true;
    problem.is_independent = true;

end


function g = disconnected_gfun(X, problem)

    x1 = X(:,1);
    x2 = X(:,2);

    a = problem.a;
    b = problem.b;
    c = problem.c;

    epsilon = 0.05;  % small asymmetry
    g1 = (x2.^2 + a*(x1 - c).^2) - b;
    g2 = (x2.^2 + a*(x1 + c).^2) - b + epsilon;
    
    k = 10; % smoothing parameter (try 5–20)
    g = -log(exp(-k*g1) + exp(-k*g2)) / k;

end