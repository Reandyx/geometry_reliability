function results = run_rare_event_study()

rng(1);

% ========================
% CONFIG
% ========================
betas = [3, 4, 5, 6];
N = 1e5;

results = struct();

for i = 1:length(betas)

beta_target = betas(i);

% --- problem ---
problem = get_problem_synthetic_rare(beta_target);

% ========================
% METHODS
% ========================
res_form = run_form(problem);
res_mcs  = run_mcs(problem, N);

% ========================
% SAFE IMPORTANCE SAMPLING
% ========================
if ~res_form.converged || isempty(res_form.U_star) || any(isnan(res_form.U_star))

    warning('IS skipped: invalid FORM design point');

    res_is = struct();
    res_is.Pf   = NaN;
    res_is.COV  = NaN;
    res_is.flag = 0;

else

    % --- FORCE CORRECT SHAPE ---
    u_star = res_form.U_star;
    u_star = double(u_star(:)');   % ensure 1×d row vector

    try
        res_is = run_is(problem, u_star, N);
    catch
        warning('IS failed: mvnrnd input issue');

        res_is = struct();
        res_is.Pf   = NaN;
        res_is.COV  = NaN;
        res_is.flag = 0;
    end

end

% ========================
% STORE RAW RESULTS ONLY
% ========================
results(i).beta = beta_target;
results(i).FORM = res_form;
results(i).MCS  = res_mcs;
results(i).IS   = res_is;

end

% ========================
% SAVE
% ========================
if ~exist('results/rare_event/raw', 'dir')
mkdir('results/rare_event/raw');
end

save('results/rare_event/raw/data.mat', 'results');

end
