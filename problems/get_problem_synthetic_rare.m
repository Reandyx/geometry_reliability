function problem = get_problem_synthetic_rare(d, b)

    if nargin < 1
        d = 2.0;   % shifts design point → controls beta
    end
    if nargin < 2
        b = 1.0;   % threshold
    end

    problem.name = sprintf('Synthetic Rare Event (d=%.2f)', d);
    problem.model_type = 'rare';

    problem.dimension = 2;

    problem.mu = [0, 0];
    problem.sigma = [1, 1];
    problem.dist = {'normal','normal'};

    problem.d = d;
    problem.b = b;

    problem.gfun = @(X) rare_gfun(X, problem);
    problem.sample = @(N) sample_inputs(problem, N);

    problem.is_normal = true;
    problem.is_independent = true;

end


function g = rare_gfun(X, problem)

    x1 = X(:,1);
    x2 = X(:,2);

    d = problem.d;
    b = problem.b;

    % KEY: shifted + slight asymmetry (fixes gradient issue)
    g = ( (x1 - d).^2 + x2.^2 + 0.1*x2 ) - b;

end