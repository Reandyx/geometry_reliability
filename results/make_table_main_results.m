function table = make_table_main_results()

% =========================
% LOAD DATA
% =========================
S = load('results/curvature/tables/summary.mat');
curv = S.summary;

S = load('results/topology/tables/summary.mat');
topo = S.summary;

S = load('results/rare_event/tables/summary.mat');
rare = S.summary;

% =========================
% CURVATURE THRESHOLDS
% =========================
gamma_10 = curv.gamma_10;
gamma_20 = curv.gamma_20;
gamma_30 = curv.gamma_30;

% =========================
% TOPOLOGY METRIC
% =========================
topology_failure = min(topo.capture_ratio);

% =========================
% RARE EVENT LIMIT
% =========================
valid_idx = rare.N_required < 1e8;

if any(valid_idx)
    beta_limit = max(rare.beta(valid_idx));
else
    beta_limit = NaN;
end

% =========================
% BUILD TABLE
% =========================
table = struct();

table.gamma_10 = gamma_10;
table.gamma_20 = gamma_20;
table.gamma_30 = gamma_30;

table.topology_failure = topology_failure;

table.beta_limit = beta_limit;

% =========================
% SAVE
% =========================
if ~exist('results/tables','dir')
    mkdir('results/tables');
end

save('results/tables/main_table.mat','table');

% =========================
% PRINT (VERY USEFUL)
% =========================
fprintf('\n============================\n');
fprintf(' MAIN RESULTS TABLE\n');
fprintf('============================\n');

fprintf('gamma_10 (10%% error) = %.4f\n', gamma_10);
fprintf('gamma_20 (20%% error) = %.4f\n', gamma_20);
fprintf('gamma_30 (30%% error) = %.4f\n', gamma_30);
fprintf('Topology min capture = %.4f\n', topology_failure);
fprintf('Beta limit (MCS feasible) = %.2f\n', beta_limit);

end