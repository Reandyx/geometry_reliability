function summary = process_curvature_results()

load('results/curvature/raw/data.mat');

n = length(results);

% --- preallocate ---
gamma = zeros(1,n);
Pf_form = zeros(1,n);
Pf_sorm = zeros(1,n);
Pf_mcs  = zeros(1,n);
err_form = zeros(1,n);
err_sorm = zeros(1,n);

for i = 1:n

    gamma(i) = results(i).gamma;

    Pf_form(i) = results(i).FORM.Pf;
    Pf_sorm(i) = results(i).SORM.Pf;
    Pf_mcs(i)  = results(i).MCS.Pf;

    err_form(i) = abs(Pf_form(i) - Pf_mcs(i)) / Pf_mcs(i);
    err_sorm(i) = abs(Pf_sorm(i) - Pf_mcs(i)) / Pf_mcs(i);

end

% SORT BY GAMMA 
[gamma_sorted, idx] = sort(gamma);

err_form_sorted = err_form(idx);
err_sorm_sorted = err_sorm(idx);

% STORE
summary.gamma = gamma;
summary.gamma_sorted = gamma_sorted;

summary.err_form = err_form;
summary.err_sorm = err_sorm;

summary.err_form_sorted = err_form_sorted;
summary.err_sorm_sorted = err_sorm_sorted;

summary.improvement = err_form_sorted - err_sorm_sorted;

% THRESHOLDS
gamma_10 = gamma_sorted(find(err_form_sorted > 0.10, 1));
gamma_20 = gamma_sorted(find(err_form_sorted > 0.20, 1));
gamma_30 = gamma_sorted(find(err_form_sorted > 0.30, 1));

if isempty(gamma_10), gamma_10 = NaN; end
if isempty(gamma_20), gamma_20 = NaN; end
if isempty(gamma_30), gamma_30 = NaN; end

summary.gamma_10 = gamma_10;
summary.gamma_20 = gamma_20;
summary.gamma_30 = gamma_30;

% SAVE
if ~exist('results/curvature/tables', 'dir')
    mkdir('results/curvature/tables');
end

save('results/curvature/tables/summary.mat', 'summary');

end