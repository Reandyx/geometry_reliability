function plot_rare_event_cost()

load('results/rare_event/tables/summary.mat');

figure;
semilogy(summary.beta, summary.N_required, '-o');

xlabel('\beta');
ylabel('Required Samples (log scale)');
set(gca,'FontSize',12);
grid on;

end