function plot_method_feasibility_map()

load('results/final/curvature_regimes.mat');
load('results/final/rare_event_regimes.mat');

g_low    = curvature_regimes.low_max;
g_high   = curvature_regimes.high_min;
beta_lim = rare_regime.beta_limit;
beta_trans = rare_regime.beta_transition;

gamma_vals = linspace(0,0.6,200);
beta_vals  = linspace(2,6,200);

[Gamma,Beta] = meshgrid(gamma_vals,beta_vals);

k = 30;

s1 = 1 ./ (1 + exp(-k*(Gamma - g_low)));
s2 = 1 ./ (1 + exp(-k*(Gamma - g_high)));
sb = 1 ./ (1 + exp(-k*(Beta - beta_lim)));

Z = 1*(1-s1) + ...
    2*(s1.*(1-s2)) + ...
    3*(s2.*(1-sb)) + ...
    4*(sb);

figure;

contourf(Gamma, Beta, Z, 20, 'LineColor','none');
hold on;

colormap(parula);
caxis([1.2 3.8]);

% Boundaries
xline(g_low,  '-k', 'LineWidth',2.5);
xline(g_high, '-k', 'LineWidth',2.5);
yline(beta_lim, '-k', 'LineWidth',2.5);     % β = 5
yline(beta_trans, '--k', 'LineWidth',1.8);  % β = 4

xlabel('\gamma');
ylabel('\beta');
title({'Method Feasibility Map', ...
       'Based on curvature, topology, rare-event'});

% Labels
text(0.05, 3.1, 'FORM','FontWeight','bold','HorizontalAlignment','center');
text(0.2, 3.1, 'SORM','FontWeight','bold','HorizontalAlignment','center');
text(0.5, 3.1, 'SORM + validation','FontWeight','bold','HorizontalAlignment','center');

text(0.3, 4.3, 'Transition regime', ...
    'FontAngle','italic','HorizontalAlignment','center');

text(0.3, 5.5, 'Importance Sampling', ...
    'FontWeight','bold','HorizontalAlignment','center');

text(0.02, beta_trans+0.1, 'Transition (high CoV)', 'FontSize',10);

box on;
grid on;
end