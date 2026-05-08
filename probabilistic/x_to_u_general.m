function U = x_to_u_general(X, problem)

% Ensure matrix shape
X = reshape(X, [], length(problem.mu));

U = zeros(size(X));

for i = 1:length(problem.mu)

    dist = lower(problem.dist{i});

    mu = problem.mu(i);
    sigma = problem.sigma(i);

    switch dist

        % =====================================
        % NORMAL
        % =====================================
        case 'normal'

            P = normcdf( (X(:,i) - mu) ./ sigma );

        % =====================================
        % LOGNORMAL
        % =====================================
        case 'lognormal'

            sigma_ln = sqrt(log(1 + (sigma/mu)^2));

            mu_ln = log(mu) - 0.5*sigma_ln^2;

            P = logncdf(X(:,i), mu_ln, sigma_ln);

        % =====================================
        % GUMBEL
        % =====================================
        case 'gumbel'

            beta_g = sigma * sqrt(6)/pi;

            mu_g = mu - 0.5772156649 * beta_g;

            z = (X(:,i) - mu_g) ./ beta_g;

            P = exp(-exp(-z));

        otherwise

            error('Unsupported distribution');

    end

    % =====================================
    % Numerical safeguard
    % =====================================
    P = min(max(P, 1e-12), 1 - 1e-12);

    % =====================================
    % Transform to U-space
    % =====================================
    U(:,i) = norminv(P);

end

end