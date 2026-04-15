function g = evaluate_limit_state_u(problem, U)

    U = reshape(U, 1, []);   % force row
    X = u_to_x(U, problem);

    g = evaluate_limit_state(problem, X);
    g = g(1);

end