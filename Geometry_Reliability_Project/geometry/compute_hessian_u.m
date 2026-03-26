function H = compute_hessian_u(problem, U)

% COMPUTE_HESSIAN_U
% Central finite difference Hessian in U-space (scalar-safe)

dim = length(U);
H = zeros(dim);

h = 1e-4;   % Larger than gradient step (good practice)

% Base evaluation (force scalar)
g0 = evaluate_limit_state_u(problem, U);
g0 = g0(1);

for i = 1:dim
    for j = i:dim
        
        if i == j
            % --- Diagonal term ---
            U_forward = U;
            U_backward = U;
            
            U_forward(i) = U_forward(i) + h;
            U_backward(i) = U_backward(i) - h;
            
            g_f = evaluate_limit_state_u(problem, U_forward);
            g_b = evaluate_limit_state_u(problem, U_backward);
            
            % Force scalar
            g_f = g_f(1);
            g_b = g_b(1);
            
            H(i,i) = (g_f - 2*g0 + g_b) / (h^2);
            
        else
            % --- Mixed partial derivative ---
            
            U_pp = U; U_pp([i j]) = U_pp([i j]) + h;
            U_pm = U; U_pm(i) = U_pm(i) + h; U_pm(j) = U_pm(j) - h;
            U_mp = U; U_mp(i) = U_mp(i) - h; U_mp(j) = U_mp(j) + h;
            U_mm = U; U_mm([i j]) = U_mm([i j]) - h;
            
            g_pp = evaluate_limit_state_u(problem, U_pp);
            g_pm = evaluate_limit_state_u(problem, U_pm);
            g_mp = evaluate_limit_state_u(problem, U_mp);
            g_mm = evaluate_limit_state_u(problem, U_mm);
            
            % Force scalar
            g_pp = g_pp(1);
            g_pm = g_pm(1);
            g_mp = g_mp(1);
            g_mm = g_mm(1);
            
            H(i,j) = (g_pp - g_pm - g_mp + g_mm) / (4*h^2);
            H(j,i) = H(i,j); % symmetry
        end
        
    end
end
X_star = u_to_x(U, problem);

disp('U ='); disp(U)
disp('X_star ='); disp(X_star)
% Enforce symmetry (numerical stability)
H = 0.5 * (H + H');

end