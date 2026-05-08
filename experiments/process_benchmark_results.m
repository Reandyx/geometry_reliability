function summary = process_benchmark_results()

load('results/benchmark/raw/data.mat');

n = length(results);

for i = 1:n

    summary.name{i} = results(i).name;
    summary.gamma(i) = results(i).gamma;
    summary.beta(i)  = results(i).beta;

    summary.err_form(i) = results(i).err_form;
    summary.err_sorm(i) = results(i).err_sorm;

    summary.predicted{i} = results(i).predicted;
    summary.best{i} = results(i).best;

    summary.correct(i) = strcmp(results(i).predicted, results(i).best);

end

summary.accuracy = mean(summary.correct);

save('results/benchmark/tables/summary.mat','summary');

end