function plot_topology_error()

load('results/topology/tables/summary.mat');

figure;
plot(summary.c, summary.error_form, '-o');
yline(1,'--','Perfect prediction');

xlabel('c (topology parameter)');
ylabel('|P_f^{FORM} - P_f^{MCS}| / P_f^{MCS}')
set(gca,'FontSize',12);
grid on;

end