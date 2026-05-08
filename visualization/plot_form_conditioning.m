function plot_form_conditioning(results)

if nargin < 1

    data = load('results/conditioning/raw/data.mat');

    results = data.results;

end

gamma = [results.gamma];

[gamma, idx] = sort(gamma);

condH = [results(idx).cond_H];
iterations = [results(idx).iterations];

form_error = [results(idx).form_error];
sorm_error = [results(idx).sorm_error];

% =====================================
% Plot 1 — Conditioning
% =====================================
figure;

semilogy(gamma, condH, '-o', 'LineWidth', 2);

xlabel('\gamma');
ylabel('cond(H)');
title('Hessian Conditioning vs Curvature');

grid on;

% =====================================
% Plot 2 — Solver Difficulty
% =====================================
figure;

plot(gamma, iterations, '-o', 'LineWidth', 2);

xlabel('\gamma');
ylabel('FORM iterations');

title('FORM Iteration Count vs Curvature');

xline(0.1,'--k');
xline(0.4,'--k');
xline(1.0,'--k');

grid on;

% =====================================
% Plot 3 — Approximation Error
% =====================================
figure;

plot(gamma, form_error, '-o', 'LineWidth', 2);
hold on;

plot(gamma, sorm_error, '-s', 'LineWidth', 2);

xlabel('\gamma');
ylabel('Relative error');

title('Approximation Error vs Curvature');

legend('FORM','SORM');

xline(0.1,'--k');
xline(0.4,'--k');
xline(1.0,'--k');

grid on;

end