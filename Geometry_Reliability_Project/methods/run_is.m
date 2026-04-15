function [Pf_IS, CoV_IS, results] = run_is(problem, N, u_star)

% ========================================
% IMPORTANCE SAMPLING (ROBUST VERSION)
% ========================================

dim = problem.dimension;

% --- IS distribution (centered at FORM design point) ---
mu_is = u_star(:)';                 % ensure row vector

% --- IS distribution (scaled) ---
sigma_scale = 1.5;   % tune: 1.2 → 2.0

Sigma_is = (sigma_scale^2) * eye(dim);

U = mvnrnd(mu_is, Sigma_is, N);

% --- Transform to X-space (standard normal → physical) ---
X = U;   % since problem.is_normal = true (no transform needed)

% --- Evaluate limit-state ---
g = problem.gfun(X);

% --- Indicator function ---
I = (g <= 0);

% ========================================
% PDF EVALUATION (LOG-SPACE)
% ========================================

% Standard normal PDF (target)
log_pdf_target = -0.5 * sum(U.^2, 2) - (dim/2)*log(2*pi);

% IS PDF (shifted normal)
diff = U - mu_is;

log_pdf_is = -0.5 * sum(diff.^2, 2) / (sigma_scale^2) ...
             - (dim/2)*log(2*pi) ...
             - dim*log(sigma_scale);

% --- Log-weights (STABLE) ---
log_w = log_pdf_target - log_pdf_is;

% Stabilization trick (VERY IMPORTANT)
log_w = log_w - max(log_w);

w = exp(log_w);

% ========================================
% ESTIMATION
% ========================================

Pf_IS = sum(w .* I) / sum(w);

% Variance estimate
var_IS = sum(w.^2 .* (I - Pf_IS).^2) / (sum(w)^2);

CoV_IS = sqrt(var_IS) / Pf_IS;

% ========================================
% DIAGNOSTICS 
% ========================================

% Effective Sample Size (ESS)
ESS = (sum(w)^2) / sum(w.^2);
ESS_ratio = ESS / N;

% Failure counts
fail_IS = sum(I);

% Weighted failure contribution
weighted_fail = sum(w .* I);

% ========================================
% PRINT (for debugging / analysis)
% ========================================

fprintf('[IS] Pf = %.3e | CoV = %.3f | ESS = %.0f (%.1f%%) | Failures = %d\n', ...
    Pf_IS, CoV_IS, ESS, 100*ESS_ratio, fail_IS);

% ========================================
% OUTPUT STRUCT
% ========================================

results = struct();
results.Pf = Pf_IS;
results.CoV = CoV_IS;
results.ESS = ESS;
results.ESS_ratio = ESS_ratio;
results.failures = fail_IS;
results.weighted_fail = weighted_fail;
results.mean_weight = mean(w);
results.max_weight = max(w);
results.min_weight = min(w);

end