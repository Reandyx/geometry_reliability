function method = select_method(gamma, beta, capture, framework)

g_low    = framework.parameters.gamma_low;
g_high   = framework.parameters.gamma_high;
beta_lim   = framework.parameters.beta_limit;
beta_trans = framework.parameters.beta_transition;
cap_thr  = framework.parameters.capture_threshold;

% ==============================
% PRIORITY 1 — TOPOLOGY
% ==============================
if capture < cap_thr
    method = 'MCS / multi-point FORM';
    return;
end

% ==============================
% PRIORITY 2 — RARE EVENT
% ==============================
if beta >= beta_lim
    method = 'Importance Sampling';
    return;
elseif beta >= beta_trans
    method = 'MCS (transition regime - caution)';
    return;
end

% ==============================
% PRIORITY 3 — CURVATURE
% ==============================
if gamma < g_low
    method = 'FORM';
elseif gamma < g_high
    method = 'SORM';
else
    method = 'SORM + validation';
end

end