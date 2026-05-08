function plot_distribution_sensitivity(results)

if nargin < 1

    data = load('results/distribution/raw/data.mat');

    results = data.results;

end

dist_names = {results.distribution};

beta_vals = [results.beta];
gamma_vals = [results.gamma];

delta_beta = [results.delta_beta];

Pf_ratio = [results.Pf_ratio];

cond_vals = zeros(size(results));
iter_vals = zeros(size(results));

for i = 1:length(results)

    cond_vals(i) = results(i).FORM.cond_H;
    iter_vals(i) = results(i).FORM.n_iterations;

end

% =====================================
% Plot 1 — Beta shift
% =====================================
figure;

bar(beta_vals);

set(gca,'XTickLabel',dist_names);

ylabel('\beta');

title('Reliability Index vs Distribution');

grid on;

% =====================================
% Plot 2 — Gamma shift
% =====================================
figure;

bar(gamma_vals);

set(gca,'XTickLabel',dist_names);

ylabel('\gamma');

title('Geometry Metric vs Distribution');

grid on;

% =====================================
% Plot 3 — Conditioning
% =====================================
figure;

semilogy(cond_vals,'-o','LineWidth',2);

set(gca,'XTick',1:length(dist_names));
set(gca,'XTickLabel',dist_names);

ylabel('cond(H)');

title('Hessian Conditioning vs Distribution');

grid on;

% =====================================
% Plot 4 — Iterations
% =====================================
figure;

bar(iter_vals);

set(gca,'XTickLabel',dist_names);

ylabel('Iterations');

title('FORM Iteration Count vs Distribution');

grid on;

% =====================================
% Plot 5 — Delta Beta
% =====================================
figure;

bar(delta_beta);

set(gca,'XTickLabel',dist_names);

ylabel('\Delta \beta');

title('Reliability Index Shift Relative to Normal');

grid on;

% =====================================
% Plot 6 — Pf Ratio
% =====================================
figure;

bar(Pf_ratio);

set(gca,'XTickLabel',dist_names);

ylabel('Pf / Pf_{normal}');

title('Failure Probability Ratio Relative to Normal');

grid on;

end