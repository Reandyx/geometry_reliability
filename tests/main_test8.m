clc;
clear;
close all;
addpath(genpath(pwd));
rng(1);

fprintf('========================================\n');
fprintf('   MAIN TEST 8 — RARE EVENT STUDY (M8)\n');
fprintf('========================================\n\n');

%% =========================
%% RUN EXPERIMENT
%% =========================
results = run_rare_event_study();

n = length(results);

%% =========================
%% EXTRACT DATA
%% =========================
c        = [results.c];
Pf_MC    = [results.Pf_MC];
Pf_FORM  = [results.Pf_FORM];
Pf_SORM  = [results.Pf_SORM];
CoV      = [results.CoV];
err_F    = [results.err_FORM];
err_S    = [results.err_SORM];
N        = [results.N];
failures = [results.failures];
efficiency = [results.efficiency];

% Improvement
improvement = (err_F - err_S) ./ err_F;

%% =========================
%% GEOMETRY CLASSIFICATION 
%% =========================
geom_class = cell(n,1);

for i = 1:n
    topology_flag = 0; % M8 = no topology issues

    geom_class{i} = classify_geometry( ...
        results(i).N, ...
        [], ...
        topology_flag, ...
        results(i).Pf_MC ...
    );
end

%% =========================
%% FIGURE 1 — ERROR vs Pf
%% =========================
figure;
loglog(Pf_MC, err_F, '-o', 'LineWidth', 2); hold on;
loglog(Pf_MC, err_S, '-s', 'LineWidth', 2);
grid on;
xlabel('P_f (Monte Carlo)');
ylabel('Relative Error');
legend('FORM','SORM','Location','best');
title('Error vs Probability (Rare Event Regime)');

%% =========================
%% FIGURE 2 — CoV vs Pf
%% =========================
figure;
loglog(Pf_MC, CoV, '-o', 'LineWidth', 2); hold on;
grid on;

xlabel('P_f (Monte Carlo)');
ylabel('Coefficient of Variation (MC)');
title('Monte Carlo Variance Explosion');

xline(1e-2, '--k', 'Transition');
xline(1e-4, '--r', 'Rare-event');

set(gca, 'XDir','reverse');

%% =========================
%% SCALING CHECK
%% =========================
valid = ~isnan(CoV) & CoV > 0;

p = polyfit(log(Pf_MC(valid)), log(CoV(valid)), 1);

fprintf('\n===== SCALING CHECK =====\n');
fprintf('Slope (log CoV vs log Pf) = %.3f\n', p(1));
fprintf('Expected ≈ -0.5 (theoretical MC scaling)\n');

%% =========================
%% FIGURE 3 — METHOD COMPARISON
%% =========================
figure;
loglog(c, Pf_MC, '-o', 'LineWidth', 2); hold on;
loglog(c, Pf_FORM, '-s', 'LineWidth', 2);
loglog(c, Pf_SORM, '-d', 'LineWidth', 2);
grid on;
xlabel('c (rarity parameter)');
ylabel('Failure Probability');
legend('MC','FORM','SORM','Location','best');
title('Method Comparison');

%% =========================
%% FIGURE 4 — FAILURE COUNT
%% =========================
figure;
semilogx(c, failures, '-o', 'LineWidth', 2);
grid on;
xlabel('c (rarity parameter)');
ylabel('Number of failures observed');
title('Monte Carlo Breakdown (Failure Count)');

%% =========================
%% FIGURE 5 — EFFICIENCY COLLAPSE
%% =========================
figure;
loglog(Pf_MC, efficiency, '-o', 'LineWidth', 2);
grid on;
xlabel('P_f (Monte Carlo)');
ylabel('Efficiency = 1 / (CoV^2 N)');
title('Monte Carlo Efficiency Collapse');

xline(1e-2, '--k');
xline(1e-4, '--r');

set(gca, 'XDir','reverse');

%% =========================
%% FIGURE 6 — SORM IMPROVEMENT
%% =========================
figure;
semilogx(Pf_MC, improvement, '-o', 'LineWidth', 2);
grid on;
xlabel('P_f (Monte Carlo)');
ylabel('Relative Improvement');
title('SORM Improvement over FORM');

set(gca, 'XDir','reverse');

%% =========================
%% PRINT TABLE 
%% =========================
fprintf('\n================ RESULTS TABLE ================\n');
fprintf(' c   |    Pf_MC   |   Pf_FORM  |   Pf_SORM  |  CoV   | Fail | class\n');
fprintf('------------------------------------------------------------------\n');

for i = 1:n
    fprintf('%4.1f | %10.2e | %10.2e | %10.2e | %8.3e | %5.0f | %s\n', ...
        results(i).c, ...
        results(i).Pf_MC, ...
        results(i).Pf_FORM, ...
        results(i).Pf_SORM, ...
        results(i).CoV, ...
        results(i).failures, ...
        geom_class{i});
end

%% =========================
%% INTERPRETATION 
%% =========================
fprintf('\n========================================\n');
fprintf('   RARE EVENT INTERPRETATION\n');
fprintf('========================================\n');

for i = 1:n
    
    if strcmp(geom_class{i}, 'rare_event')
        regime = 'MC BREAKDOWN (use FORM/IS)';
        
    elseif err_F(i) < 0.1
        regime = 'FORM VALID';
        
    elseif abs(err_F(i) - err_S(i)) < 1e-6
        regime = 'NO CURVATURE EFFECT (SORM useless)';
        
    elseif improvement(i) > 0.1
        regime = 'SORM USEFUL';
        
    else
        regime = 'MODEL BIAS (FORM underestimates)';
    end
    
    fprintf('c=%4.1f | Pf=%.1e | class=%s | %s\n', ...
        c(i), Pf_MC(i), geom_class{i}, regime);
end

fprintf('\n=== M8 COMPLETE ===\n');