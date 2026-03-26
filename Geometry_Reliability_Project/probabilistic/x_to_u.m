function U = x_to_u(X, problem)

    % Ensure column vector
    X = X(:);

    d = length(problem.mu);

    % Safety check
    if length(X) ~= d
        error('x_to_u: dimension mismatch');
    end

    U = zeros(d,1);

    for i = 1:d
        U(i) = (X(i) - problem.mu(i)) / problem.sigma(i);
    end

end