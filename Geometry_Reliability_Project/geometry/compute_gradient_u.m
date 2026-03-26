function grad = compute_gradient_u(problem, U)

dim = length(U);
grad = zeros(dim,1);

h = 1e-6;

for i = 1:dim
    
    U_forward = U;
    U_backward = U;

    U_forward(i) = U_forward(i) + h;
    U_backward(i) = U_backward(i) - h;

    g_forward = evaluate_limit_state_u(problem, U_forward);
    g_backward = evaluate_limit_state_u(problem, U_backward);

    % Ensure scalar (FORM assumption)
    g_forward = g_forward(1);
    g_backward = g_backward(1);

    % Stability check
    if any(~isfinite(g_forward)) || any(~isfinite(g_backward))
        warning('Gradient unstable at index %d', i);
        grad(i) = 0;
        continue;
    end

    grad(i) = (g_forward - g_backward) / (2*h);
end

end