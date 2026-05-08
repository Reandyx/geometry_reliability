function problem = get_problem_disconnected_c(b, c)

    a = 1.5;

    problem.name = sprintf('Disconnected (c=%.2f)', c);
    problem.dimension = 2;

    problem.mu = [0, 0];
    problem.sigma = [1, 1];

    problem.dist = {'normal','normal'};

    problem.a = a;
    problem.b = b;
    problem.c = c;

    % SOFTMIN SHARPNESS PARAMETER
    problem.k = 20;

    problem.gfun = @(U) disconnected_gfun_c(U, problem);
    problem.sample = @(N) sample_inputs(problem, N);

    problem.is_normal = true;
    problem.is_independent = true;

end

function g = disconnected_gfun_c(U, problem)

    x1 = U(:,1);
    x2 = U(:,2);

    a = problem.a;
    b = problem.b;
    c = problem.c;
    k = problem.k;

    epsilon = 0.05;

    % --- branch functions ---
    g1 = (x2.^2 + a*(x1 - c).^2) - b;
    g2 = (x2.^2 + a*(x1 + c).^2) - b + epsilon;

    % --- numerically stable softmin ---
    m = min([g1, g2], [], 2);

    g = m - (1/k)*log( exp(-k*(g1 - m)) + exp(-k*(g2 - m)) );

end