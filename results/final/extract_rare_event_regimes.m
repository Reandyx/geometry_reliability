function rare_regime = extract_rare_event_regimes()

load('results/rare_event/tables/summary.mat');

beta = summary.beta;
cov  = summary.CoV_mcs;
Nreq = summary.N_required;

fprintf('\n[RARE EVENT REGIMES]\n');

% =========================================
% DEFINE CONDITIONS
% =========================================

% Transition (warning): high variance or high cost
idx_transition = (cov > 0.1) | (Nreq > 1e8);

% Failure (true breakdown): solver collapses
idx_failure = isnan(cov) | isnan(Nreq);

% =========================================
% EXTRACT THRESHOLDS
% =========================================

% Transition threshold (first sign of instability)
if any(idx_transition)
    beta_transition = min(beta(idx_transition));
else
    beta_transition = NaN;
    warning('No transition regime detected');
end

% Failure threshold (true breakdown)
if any(idx_failure)
    beta_fail = min(beta(idx_failure));
else
    beta_fail = NaN;
    warning('No failure regime detected');
end

% =========================================
% SAFE REGION
% =========================================

idx_safe = ~(idx_transition | idx_failure);

if any(idx_safe)
    beta_safe_max = max(beta(idx_safe));
else
    beta_safe_max = NaN;
end

% =========================================
% STORE RESULTS
% =========================================

rare_regime.beta_transition = beta_transition;
rare_regime.beta_limit      = beta_fail;
rare_regime.max_beta_mcs    = beta_safe_max;

save('results/final/rare_event_regimes.mat','rare_regime');

% =========================================
% PRINT
% =========================================

fprintf('Beta transition (warning) = %.2f\n', beta_transition);
fprintf('Beta limit (failure)     = %.2f\n', beta_fail);
fprintf('Max beta safe (MCS)      = %.2f\n', beta_safe_max);

end