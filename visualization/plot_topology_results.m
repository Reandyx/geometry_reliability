function plot_topology_results()

load('results/topology/tables/summary.mat');

figure;
plot(summary.c, summary.capture_ratio, '-o');
yline(1,'--','Perfect prediction');

xlabel('c');
ylabel('Capture Ratio (FORM / MCS)');
set(gca,'FontSize',12);
grid on;

end