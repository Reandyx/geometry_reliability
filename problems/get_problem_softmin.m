function problem = get_problem_softmin(beta1, beta2, c, k)

    problem.dimension = 2;

    problem.dist = {'normal','normal'};

    problem.mu = [0; 0];
    problem.sigma = [1; 1];

    problem.g = @(U) g_soft(U, beta1, beta2, c, k);
    problem.gfun = @(U) g_soft(U, beta1, beta2, c, k);

    problem.grad_g = @(U) numerical_grad(@(X) g_soft(X, beta1, beta2, c, k), U);
    problem.hess_g = @(U) numerical_hess(@(X) g_soft(X, beta1, beta2, c, k), U);

end


% ============================================================
% SOFT-MIN LIMIT STATE
% ============================================================
function val = g_soft(U, beta1, beta2, c, k)

    if isvector(U)
        U = U(:)';
    end

    g1 = beta1 - sqrt((U(:,1)-c).^2 + U(:,2).^2);
    g2 = beta2 - sqrt((U(:,1)+c).^2 + U(:,2).^2);

    m = min(g1, g2);

    val = m - (1/k)*log( exp(-k*(g1 - m)) + exp(-k*(g2 - m)) );

end

% ============================================================
% NUMERICAL GRADIENT
% ============================================================
function grad = numerical_grad(fun, U)

    % Ensure row vector
    U = U(:)';  

    n = length(U);
    grad = zeros(1,n);

    for i = 1:n

        % Adaptive step 
        h = 1e-5 * max(1, abs(U(i)));

        U_f = U;
        U_b = U;

        U_f(i) = U_f(i) + h;
        U_b(i) = U_b(i) - h;

        g_f = fun(U_f);
        g_b = fun(U_b);

        grad(i) = (g_f - g_b) / (2*h);
    end

end


% ============================================================
% NUMERICAL HESSIAN
% ============================================================
function H = numerical_hess(fun, U)

    eps = 1e-4;
    n = length(U);
    H = zeros(n);

    for i = 1:n
        for j = 1:n

            U1 = U; U2 = U; U3 = U; U4 = U;

            U1(i) = U1(i) + eps; U1(j) = U1(j) + eps;
            U2(i) = U2(i) + eps; U2(j) = U2(j) - eps;
            U3(i) = U3(i) - eps; U3(j) = U3(j) + eps;
            U4(i) = U4(i) - eps; U4(j) = U4(j) - eps;

            H(i,j) = (fun(U1) - fun(U2) - fun(U3) + fun(U4)) / (4*eps^2);

        end
    end

end