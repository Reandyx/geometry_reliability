function framework = build_method_selection_framework()

% =========================================
% LOAD REGIMES
% =========================================
load('results/final/curvature_regimes.mat');
load('results/final/topology_regimes.mat');
load('results/final/rare_event_regimes.mat');

g_low    = curvature_regimes.low_max;
g_high   = curvature_regimes.high_min;
beta_lim = rare_regime.beta_limit;

% =========================================
% BUILD RULE TABLE (EXPLICIT)
% =========================================
framework.rules = {

    sprintf('gamma < %.2f & connected', g_low),          'FORM';
    sprintf('%.2f <= gamma < %.2f', g_low, g_high),      'SORM';
    sprintf('gamma >= %.2f', g_high),                    'SORM + validation';
    sprintf('capture < %.2f (topology failure)', 0.8),   'MCS / multi-point FORM';
    sprintf('beta >= %.1f', beta_lim),                   'Importance Sampling'
    sprintf('%.1f <= beta < %.1f', beta_trans, beta_lim), 'MCS (caution)'

};

% =========================================
% STORE PARAMETERS
% =========================================
framework.parameters.gamma_low  = g_low;
framework.parameters.gamma_high = g_high;
framework.parameters.beta_limit = beta_lim;
framework.parameters.beta_transition = rare_regime.beta_transition;
framework.parameters.capture_threshold = 0.8;

framework.description = 'Geometry-aware method selection based on curvature (gamma), topology (capture), and rare-event regime (beta).';

% =========================================
% SAVE
% =========================================
save('results/final/method_framework.mat','framework');

disp('Method selection framework built (final).');

end
