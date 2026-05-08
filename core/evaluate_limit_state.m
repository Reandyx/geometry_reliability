function g = evaluate_limit_state(problem, X)

    % Force row format for single sample
    if isvector(X)
        X = reshape(X, 1, []);
    end

    if size(X,2) ~= problem.dimension
        error('Input dimension mismatch');
    end

    eval_counter('add', size(X,1));

    g = problem.gfun(X);

end