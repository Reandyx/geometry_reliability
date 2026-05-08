function plot_gamma_thresholds(results)

if nargin < 1

    data = load('results/gamma_thresholds/raw/data.mat');

    results = data.results;

end

if ~exist('results/gamma_thresholds/figures','dir')

    mkdir('results/gamma_thresholds/figures');

end

gamma = [results.gamma];

form_error = [results.form_error];

sorm_error = [results.sorm_error];

sorm_gain = [results.sorm_gain];

% =====================================
% Plot 1 — FORM error vs gamma
% =====================================
figure;

scatter(gamma, form_error, ...
    80, 'filled');

hold on;

xline(0.1,'g--','LineWidth',2);

xline(1.0,'r--','LineWidth',2);

xlabel('\gamma');

ylabel('Relative FORM error');

title('FORM Error vs Geometry Index');

grid on;

saveas(gcf, ...
'results/gamma_thresholds/figures/form_error_vs_gamma.png');

% =====================================
% Plot 2 — SORM error vs gamma
% =====================================
figure;

scatter(gamma, sorm_error, ...
    80, 'filled');

hold on;

xline(0.1,'g--','LineWidth',2);

xline(1.0,'r--','LineWidth',2);

xlabel('\gamma');

ylabel('Relative SORM error');

title('SORM Error vs Geometry Index');

grid on;

saveas(gcf, ...
'results/gamma_thresholds/figures/sorm_error_vs_gamma.png');

% =====================================
% Plot 3 — SORM gain
% =====================================
figure;

scatter(gamma, sorm_gain, ...
    80, 'filled');

hold on;

yline(0,'k--');

xlabel('\gamma');

ylabel('SORM gain');

title('SORM Improvement vs Geometry Index');

grid on;

saveas(gcf, ...
'results/gamma_thresholds/figures/sorm_gain_vs_gamma.png');

end