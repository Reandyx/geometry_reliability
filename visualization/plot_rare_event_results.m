function plot_rare_event_results()

load('results/rare_event/tables/summary.mat');

figure;
plot(summary.beta, summary.CoV_mcs, '-o');
hold on;
plot(summary.beta, summary.CoV_is, '-x');

xlabel('\beta');
ylabel('CoV');
legend('MCS', 'IS');
set(gca,'FontSize',12);
grid on;

end