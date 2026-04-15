function [kappa_global, N] = compute_curvature_metrics(H, grad_g)
% COMPUTE_CURVATURE_METRICS
% Consistent curvature metrics based on principal curvatures

    % Get principal curvatures
    kappa = compute_principal_curvatures(H, grad_g);

    if isempty(kappa)
        kappa_global = 0;
    else
        % RMS curvature (physically meaningful scalar)
        kappa_global = norm(kappa);
    end

    grad_norm = norm(grad_g);

    if grad_norm < 1e-12
        N = NaN;
    else
        % Normalized curvature index
        N = kappa_global;
    end

end