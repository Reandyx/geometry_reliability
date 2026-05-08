function X = u_to_x_general(U, problem)

U = reshape(U, [], length(problem.mu));

X = zeros(size(U));

P = normcdf(U);

for i = 1:length(problem.mu)

    dist = lower(problem.dist{i});

    mu = problem.mu(i);
    sigma = problem.sigma(i);

    switch dist

        case 'normal'

            X(:,i) = mu + sigma .* U(:,i);

        case 'lognormal'

            % convert mean/std -> lognormal params
            sigma_ln = sqrt(log(1 + (sigma/mu)^2));

            mu_ln = log(mu) - 0.5*sigma_ln^2;

            X(:,i) = logninv(P(:,i), mu_ln, sigma_ln);

        case 'gumbel'

            beta_g = sigma * sqrt(6)/pi;

            mu_g = mu - 0.5772156649 * beta_g;

            X(:,i) = mu_g - beta_g .* log(-log(P(:,i)));

        otherwise

            error('Unsupported distribution');

    end

end

end