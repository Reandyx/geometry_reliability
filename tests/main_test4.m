clc;
clear;
close all;
addpath(genpath(pwd));
rng(1);

fprintf('\n========================================\n');
fprintf('   MAIN TEST 4 — CURVATURE STUDY (M4)\n');
fprintf('========================================\n');

results_global = run_curvature_study('global');
results_local  = run_curvature_study('local');

T_global = struct2table(results_global);
T_local  = struct2table(results_local);

fprintf('\n============================\n');
fprintf('   GLOBAL RESULTS\n');
fprintf('============================\n');
disp(T_global);

fprintf('\n============================\n');
fprintf('   LOCAL RESULTS\n');
fprintf('============================\n');
disp(T_local);

figure; hold on; grid on;
plot([results_global.a], [results_global.error_FORM], '-o', 'LineWidth', 2);
plot([results_local.a],  [results_local.error_FORM],  '-s', 'LineWidth', 2);
xlabel('a');
ylabel('Relative Error');
legend('Global', 'Local');
title('FORM Error vs Curvature');

figure; hold on; grid on;
plot([results_global.a], [results_global.bias], '-o', 'LineWidth', 2);
plot([results_local.a],  [results_local.bias],  '-s', 'LineWidth', 2);
xlabel('a');
ylabel('Bias');
legend('Global', 'Local');
title('FORM Bias vs Curvature');

figure; hold on; grid on;
plot([results_global.a], [results_global.beta_FORM], '-o', 'LineWidth', 2);
plot([results_local.a],  [results_local.beta_FORM],  '-s', 'LineWidth', 2);
xlabel('a');
ylabel('\beta');
legend('Global', 'Local');
title('Reliability Index vs Curvature');

figure; hold on; grid on;
errorbar([results_global.a], [results_global.Pf_MC], [results_global.sigma_MC], '-o', 'LineWidth', 2);
plot([results_global.a], [results_global.Pf_FORM], '-o', 'LineWidth', 2);

errorbar([results_local.a], [results_local.Pf_MC], [results_local.sigma_MC], '-s', 'LineWidth', 2);
plot([results_local.a], [results_local.Pf_FORM], '-s', 'LineWidth', 2);

xlabel('a');
ylabel('Failure Probability');
legend('MC Global','FORM Global','MC Local','FORM Local');
title('Pf Comparison');

figure; hold on; grid on;
U_global = reshape([results_global.U_star], 2, []);
U_local  = reshape([results_local.U_star],  2, []);
plot(U_global(1,:), U_global(2,:), '-o', 'LineWidth', 2);
plot(U_local(1,:),  U_local(2,:),  '-s', 'LineWidth', 2);
xlabel('U_1^*');
ylabel('U_2^*');
legend('Global','Local');
title('Design Point Trajectory');

fprintf('\n========================================\n');
fprintf('   M4 COMPLETE\n');
fprintf('========================================\n');