function problem = get_problem_nonlinear()

problem.name = 'Nonlinear response model';
problem.model_type = 'engineering';

problem.dimension = 2;

problem.mu = [0, 0];
problem.sigma = [1, 1];
problem.dist = {'normal','normal'};

problem.gfun = @(X) nonlinear_gfun(X);
problem.sample = @(N) sample_inputs(problem, N);

problem.is_normal = true;
problem.is_independent = true;

end


function g = nonlinear_gfun(X)

x1 = X(:,1);
x2 = X(:,2);

g = 3 - (x1.^2 + 0.5*x2.^2 + 0.5*sin(2*pi*x1));

end