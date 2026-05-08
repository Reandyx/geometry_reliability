function plot_curvature_results_unsorted()

load('results/curvature/tables/summary.mat');

figure;
plot(summary.gamma, summary.err_form, 'o-','LineWidth',1.5);
hold on;
plot(summary.gamma, summary.err_sorm, 'x-','LineWidth',1.5);
ylim([0 0.5])

xlabel('\gamma');
ylabel('Relative Error');
legend('FORM', 'SORM');
title('Unsorted (non-physical ordering)')
set(gca,'FontSize',12);
grid on;

end