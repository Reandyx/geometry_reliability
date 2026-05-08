function topology_regime = extract_topology_regimes()

load('results/topology/tables/summary.mat');

capture = summary.capture_ratio;
error   = summary.error_form;
c       = summary.c;

fprintf('\n[TOPOLOGY REGIMES]\n');

% Failure condition
idx_fail = capture < 0.8;

if any(idx_fail)
    c_fail = min(c(idx_fail));
else
    c_fail = NaN;
    warning('No topology failure detected');
end

topology_regime.failure_threshold = c_fail;
topology_regime.min_capture = min(capture);

save('results/final/topology_regimes.mat','topology_regime');

fprintf('Failure threshold c = %.4f\n', c_fail);
fprintf('Minimum capture = %.4f\n', min(capture));

end