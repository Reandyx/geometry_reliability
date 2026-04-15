function X = sample_inputs(problem, N)

    d = problem.dimension;
    X = zeros(N, d);

    for i = 1:d
        if strcmp(problem.dist{i}, 'normal')
            X(:,i) = problem.mu(i) + problem.sigma(i) .* randn(N,1);
        else
            error('Unsupported distribution');
        end
    end

end