function summary = process_topology_results()

load('results/topology/raw/data.mat');

n = length(results);

for i = 1:n

c(i) = results(i).c;

Pf_form(i) = results(i).FORM.Pf;
Pf_mcs(i)  = results(i).MCS.Pf;

capture_ratio(i) = results(i).capture_ratio;
discrepancy(i) = abs(Pf_form(i) - Pf_mcs(i));
err_form(i) = abs(Pf_form(i) - Pf_mcs(i)) / Pf_mcs(i);

end

summary.c = c;
summary.capture_ratio = capture_ratio;
summary.discrepancy = discrepancy;
summary.error_form = err_form;

if ~exist('results/topology/tables', 'dir')
mkdir('results/topology/tables');
end

save('results/topology/tables/summary.mat', 'summary');

end
