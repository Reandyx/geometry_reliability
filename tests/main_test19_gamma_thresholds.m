clc;
clear;
close all;

fprintf('\n');
fprintf('=====================================\n');
fprintf('MAIN TEST 16 — GAMMA THRESHOLDS\n');
fprintf('=====================================\n');

% =====================================
% Run study
% =====================================
results = run_gamma_threshold_study();

% =====================================
% Generate plots
% =====================================
plot_gamma_thresholds(results);

% =====================================
% Summary
% =====================================
fprintf('\n');
fprintf('=====================================\n');
fprintf('GAMMA THRESHOLD SUMMARY\n');
fprintf('=====================================\n');

for i = 1:length(results)

    fprintf('\n');

    fprintf('a = %6.2f\n', ...
        results(i).a);

    fprintf('gamma = %8.4f\n', ...
        results(i).gamma);

    fprintf('FORM err = %8.4f\n', ...
        results(i).form_error);

    fprintf('SORM err = %8.4f\n', ...
        results(i).sorm_error);

    fprintf('SORM gain = %8.4f\n', ...
        results(i).sorm_gain);

    fprintf('regime = %s\n', ...
        results(i).regime);

end

fprintf('\n');
fprintf('MAIN TEST 19 COMPLETED.\n');
fprintf('\n');