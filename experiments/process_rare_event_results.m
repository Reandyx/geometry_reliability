function summary = process_rare_event_results()

load('results/rare_event/raw/data.mat');

n = length(results);

for i = 1:n
    
    beta(i) = results(i).beta;
    
    Pf_mcs(i) = results(i).MCS.Pf;
    
    if isfield(results(i).MCS, 'COV')
        CoV_mcs(i) = results(i).MCS.COV;
    else
        CoV_mcs(i) = NaN;
    end
    
    if isfield(results(i).IS, 'COV')
        CoV_is(i) = results(i).IS.COV;
    else
        CoV_is(i) = NaN;
    end

    Pf = results(i).MCS.Pf;
    
    if Pf > 0
        N_required(i) = 100 / Pf;
    else
        N_required(i) = NaN;
    end

end

summary.beta = beta;
summary.CoV_mcs = CoV_mcs;
summary.CoV_is = CoV_is;
summary.N_required = N_required;

if ~exist('results/rare_event/tables', 'dir')
mkdir('results/rare_event/tables');
end

save('results/rare_event/tables/summary.mat', 'summary');

end
