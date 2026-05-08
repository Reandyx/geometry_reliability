function plot_benchmark_comparison()

load('results/benchmark/tables/summary.mat');

figure;

gamma_plot = summary.gamma;
gamma_plot(gamma_plot < 1e-6) = 1e-6;

bar(gamma_plot);
hold on;

% threshold
yline(0.1, '--', '\gamma = 0.1 threshold');

set(gca,'XTickLabel', summary.name);
set(gca,'YScale','log');

ylabel('\gamma (log scale)');
title('Curvature of Benchmark Problems');

grid on;

% === ADD TEXT LABELS ===
for i = 1:length(summary.gamma)
    text(i, gamma_plot(i)*1.5, sprintf('%.2e', summary.gamma(i)), ...
        'HorizontalAlignment','center');
end

end