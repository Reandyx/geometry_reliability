function kappa = compute_principal_curvatures(H, grad_g)

    % Normalize gradient → normal direction
    alpha = grad_g / norm(grad_g);

    n = length(alpha);

    % Projection matrix onto tangent space
    P = eye(n) - alpha * alpha';

    % Projected Hessian
    H_t = P * H * P;

    [vecs, vals] = eig(H_t);

disp('--- CURVATURE DEBUG ---');
disp('H_t ='); disp(H_t);
disp('Eigenvalues ='); disp(diag(vals));
disp('Eigenvectors ='); disp(vecs);
disp('alpha ='); disp(alpha);
    
    lambda = diag(vals);
    
    % Remove eigenvalue aligned with alpha
    keep = true(length(lambda),1);
    
    for i = 1:length(lambda)
        if abs(vecs(:,i)' * alpha) > 1 - 1e-6
            keep(i) = false;
        end
    end
    
    lambda = lambda(keep);

    % Principal curvatures
    kappa = lambda(:) / norm(grad_g);

end