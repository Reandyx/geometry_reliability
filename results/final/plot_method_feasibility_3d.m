function plot_method_feasibility_3d()

load('results/final/curvature_regimes.mat');
load('results/final/rare_event_regimes.mat');

g_low    = curvature_regimes.low_max;
g_high   = curvature_regimes.high_min;
beta_lim = rare_regime.beta_limit;
beta_trans = rare_regime.beta_transition;

gamma_vals = linspace(0,0.6,150);
beta_vals  = linspace(2,6,150);

[Gamma,Beta] = meshgrid(gamma_vals,beta_vals);

% =========================================
% CONTINUOUS TRANSITIONS (SIGMOID)
% =========================================
k = 10; % sharpness

% curvature transitions
s1 = 1 ./ (1 + exp(-k*(Gamma - g_low)));
s2 = 1 ./ (1 + exp(-k*(Gamma - g_high)));

% beta transition
sb = 1 ./ (1 + exp(-k*(Beta - beta_lim)));

% =========================================
% BUILD SMOOTH REGION FIELD
% =========================================
Z = 1*(1-s1) + ...
    2*(s1.*(1-s2)) + ...
    3*(s2.*(1-sb)) + ...
    4*(sb);

% =========================================
% PLOT
% =========================================
figure;

surf(Gamma, Beta, Z, ...
    'EdgeColor','none', ...
    'FaceAlpha',0.95);

hold on;   

colormap(parula);
shading interp;

% LIGHTING
lighting phong;
camlight headlight;

% =========================
% LABELS
% =========================
xlabel('\gamma');
ylabel('\beta');
zlabel('Method index (smoothed)');

title({'Smoothed Method Selection Landscape', ...
       '(Visualization of discrete decision framework)'});

colorbar;

% =========================
% CONTOURS (clean)
% =========================
contour3(Gamma, Beta, Z, 8, ...
    'k', ...
    'LineWidth', 0.6);
end