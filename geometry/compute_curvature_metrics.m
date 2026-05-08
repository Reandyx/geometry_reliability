function metrics = compute_curvature_metrics(beta, kappa)
% COMPUTE_CURVATURE_METRICS
% Canonical curvature metric for Paper 2
%
% gamma = beta * max |kappa_i|

    if isempty(kappa)
        kappa_max = 0;
    else
        kappa_max = max(abs(kappa));
    end

    gamma = beta * kappa_max;

    metrics.kappa = kappa;
    metrics.kappa_max = kappa_max;
    metrics.gamma = gamma;

end