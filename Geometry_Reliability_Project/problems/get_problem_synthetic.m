function problem = get_problem_synthetic(a, b, c)

    % Defaults
    if nargin < 1
        a = 1.0;
    end

    if nargin < 2
        b = 4.0; % controls Pf
    end

    if nargin < 3
        c = 0.0; % controls topology / nonlinearity
    end

    % Problem definition
    problem.name = sprintf('Synthetic Geometry (a=%.2f, c=%.2f)', a, c);
    problem.model_type = 'global';
    problem.dimension = 2;

    problem.mu = [0, 0];
    problem.sigma = [1, 1];

    problem.dist = {'normal', 'normal'};

    % Parameters
    problem.a = a;
    problem.b = b;
    problem.c = c;

    % Function handles
    problem.gfun = @(X) synthetic_gfun(X, problem);
    problem.sample = @(N) sample_inputs(problem, N);

    % Flags
    problem.is_normal = true;
    problem.is_independent = true;

end

% ============================================================
% LIMIT STATE FUNCTION
% ============================================================
function g = synthetic_gfun(X, problem)

    x1 = X(:,1);
    x2 = X(:,2);

    a = problem.a;
    b = problem.b;
    c = problem.c;

    % --- base quadratic (curvature control) ---
    quad = (x1 - 1).^2 + a * (x2 - 0.5).^2;

    % --- topology term ---
    topo = h(x1, x2);

    % --- final limit state ---
    g = b - quad + c * topo;

end

% ============================================================
% TOPOLOGY FUNCTION (controlled nonlinearity)
% ============================================================
function val = h(x1, x2)

    % Smooth, controlled non-convexity
    val = sin(2*pi*x1);

end