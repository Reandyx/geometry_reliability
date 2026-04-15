function kappa = compute_principal_curvatures(H, grad_g)

    % Normalize gradient → normal direction
    grad = grad_g(:);
    alpha = grad / norm(grad);

    n = length(alpha);

    % Projection matrix onto tangent space
    P = eye(n) - alpha * alpha';

    % Projected Hessian
    H_t = P * H * P;

    % Eigenvalues
    lambda = eig(H_t);

    % Remove near-zero eigenvalue (normal direction)
    tol = 1e-6;
    kappa = lambda(abs(lambda) > tol);

    % Normalize by gradient magnitude
    kappa = kappa / norm(grad);

end