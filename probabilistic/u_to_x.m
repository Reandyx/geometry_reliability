function X = u_to_x(U, problem)

    % U: (N x d)
    U = reshape(U, [], length(problem.mu));

    X = zeros(size(U));

    for i = 1:length(problem.mu)
        X(:,i) = problem.mu(i) + problem.sigma(i) .* U(:,i);
    end

end