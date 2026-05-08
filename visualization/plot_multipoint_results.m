function plot_multipoint_results(results)

if nargin < 1

    data = load('results/multipoint/raw/data.mat');

    results = data.results;

end

if ~exist('results/multipoint/figures','dir')
    mkdir('results/multipoint/figures');
end

% =====================================
% Plot 1 — Capture comparison
% =====================================
c_vals = [results.c];

single_capture = [results.capture_single];

multi_capture = [results.capture_multi];

figure;

plot(c_vals, single_capture, '-o', 'LineWidth', 2);

hold on;

plot(c_vals, multi_capture, '-s', 'LineWidth', 2);

yline(1.0,'k--');

xlabel('Topology parameter c');

ylabel('Probability capture ratio');

title('Single vs Multi-point FORM');

legend('Single FORM', ...
       'Multi-point FORM', ...
       'Perfect capture');

grid on;

saveas(gcf, ...
'results/multipoint/figures/capture_comparison.png');

% =====================================
% Plot 2 — Branch count
% =====================================
branches = zeros(size(results));

for i = 1:length(results)

    branches(i) = results(i).MULTI.n_branches;

end

figure;

bar(c_vals, branches);

xlabel('Topology parameter c');

ylabel('Detected branches');

title('Detected MPP Branches');

grid on;

saveas(gcf, ...
'results/multipoint/figures/branch_count.png');

% =====================================
% Plot 3 — Branch contribution chart
% =====================================
for i = 1:length(results)

    branches_i = ...
        results(i).MULTI.branch_results;

    Pf_i = zeros(length(branches_i),1);

    for j = 1:length(branches_i)

        Pf_i(j) = branches_i(j).Pf;

    end

    figure;

    bar(Pf_i);

    xlabel('Branch');

    ylabel('Branch Pf contribution');

    title(sprintf( ...
        'Branch Contributions (c=%.2f)', ...
        results(i).c));

    grid on;

    saveas(gcf, sprintf( ...
        'results/multipoint/figures/branch_contributions_c_%0.2f.png', ...
        results(i).c));

end

end