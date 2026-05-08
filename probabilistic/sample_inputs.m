function X = sample_inputs(problem, N)

d = problem.dimension;

X = zeros(N, d);

for i = 1:d

    dist = lower(problem.dist{i});

    mu = problem.mu(i);
    sigma = problem.sigma(i);

    switch dist

        % =====================================
        % NORMAL
        % =====================================
        case 'normal'

            X(:,i) = mu + sigma .* randn(N,1);

        % =====================================
        % LOGNORMAL
        % =====================================
        case 'lognormal'

            sigma_ln = sqrt(log(1 + (sigma/mu)^2));

            mu_ln = log(mu) - 0.5*sigma_ln^2;

            X(:,i) = lognrnd(mu_ln, sigma_ln, N, 1);

        % =====================================
        % GUMBEL
        % =====================================
        case 'gumbel'

            beta_g = sigma * sqrt(6)/pi;

            mu_g = mu - 0.5772156649 * beta_g;

            U = rand(N,1);

            X(:,i) = mu_g - beta_g .* log(-log(U));

        otherwise

            error('Unsupported distribution');

    end

end

end