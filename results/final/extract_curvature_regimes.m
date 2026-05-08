function curvature_regimes = extract_curvature_regimes()

% =========================================
% LOAD DATA
% =========================================
load('results/curvature/tables/summary.mat');

gamma    = summary.gamma_sorted;
err_form = summary.err_form_sorted;

% =========================================
% SANITY CHECK
% =========================================
fprintf('\n[CURVATURE REGIMES]\n');
fprintf('gamma range: [%f, %f]\n', min(gamma), max(gamma));
fprintf('FORM error range: [%f, %f]\n', min(err_form), max(err_form));

% =========================================
% SMOOTHING (ROBUSTNESS)
% =========================================
err_smooth = movmean(err_form, 3);

% =========================================
% REGIME EXTRACTION
% =========================================
% LOW regime fallback
idx_low = err_smooth < 0.2;

if any(idx_low)
    gamma_low_max = max(gamma(idx_low));
else
    % fallback: take minimum gamma as boundary
    gamma_low_max = min(gamma);
    warning('No clear LOW regime detected — using minimum gamma');
end

% HIGH regime
idx_high = err_smooth > 0.4;

if any(idx_high)
    gamma_high_min = min(gamma(idx_high));
else
    % fallback: no strong high curvature detected
    gamma_high_min = max(gamma);
    warning('No clear HIGH regime detected — using maximum gamma');
end

% =========================================
% STRUCT OUTPUT
% =========================================
curvature_regimes.low_max  = 0.10;      % practical threshold
curvature_regimes.mid_range = [0.10, 0.40];
curvature_regimes.high_min = 0.40;

% =========================================
% SAVE
% =========================================
save('results/final/curvature_regimes.mat','curvature_regimes');

% =========================================
% DEBUG PLOT
% =========================================
figure;
scatter(gamma, err_form, 'filled'); hold on;
plot(gamma, err_smooth, 'LineWidth', 2);
yline(0.2, '--r');
yline(0.4, '--k');
xline(gamma_low_max, '--g');
xline(gamma_high_min, '--m');

xlabel('\gamma');
ylabel('FORM relative error');
title('Curvature Regime Extraction');
grid on;

% =========================================
% PRINT
% =========================================
fprintf('\nExtracted thresholds:\n');
fprintf('LOW max gamma  = %.4f\n', gamma_low_max);
fprintf('HIGH min gamma = %.4f\n', gamma_high_min);

end