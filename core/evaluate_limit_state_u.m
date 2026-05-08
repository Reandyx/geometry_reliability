function g = evaluate_limit_state_u(problem, U)

    U = reshape(U, 1, []);

    % =====================================
    % Use fast affine transform if all
    % marginals are Gaussian
    % =====================================
    if all(strcmpi(problem.dist, 'normal'))

        X = u_to_x(U, problem);

    else

        X = u_to_x_general(U, problem);

    end

    g = evaluate_limit_state(problem, X);

    g = g(1);

end