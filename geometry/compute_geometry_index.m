function gamma = compute_geometry_index(beta, kappa)

if isempty(kappa)

    gamma = 0;

    return;

end

gamma = beta * max(abs(kappa));

end