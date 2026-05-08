function plot_curvature_gamma()

load('results/curvature/tables/summary.mat');

figure;
plot(summary.gamma_sorted, summary.err_form_sorted, '-o');
hold on;
plot(summary.gamma_sorted, summary.err_sorm_sorted, '-x');

xlabel('\gamma = \beta max|\kappa_i|');
ylabel('Relative Error');
legend('FORM','SORM');
set(gca,'FontSize',12);
grid on;

end