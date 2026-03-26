function results = run_rare_event_study()

%% Fixed geometry
a = 1.0;

%% Control parameter (rarity increases with c)
c_values = [2,4,6,8,10,12,14,16,18,20,22,24,26,28,30];

%% Initialize results
results = struct([]);

for i = 1:length(c_values)

    c = c_values(i);

    fprintf('--- CASE c = %.2f ---\n', c);

    %% Build problem
    problem = get_problem_local_curvature(a, c);

    %% Monte Carlo scaling
    if c <= 2
        N = 1e5;
    elseif c <= 3
        N = 1e6;
    else
        N = 1e7;
    end

    %% Run methods
    res_mc   = run_mcs(problem, N);
    res_form = run_form(problem);
    res_sorm = run_sorm(problem);

    Pf_MC   = res_mc.Pf;
    Pf_FORM = res_form.Pf;
    Pf_SORM = res_sorm.Pf;

    %% Monte Carlo diagnostics
    if Pf_MC == 0
        CoV = NaN;
        efficiency = 0;
        err_FORM = NaN;
        err_SORM = NaN;
        failures = 0;
    else
        CoV = sqrt((1 - Pf_MC) / (N * Pf_MC));
        efficiency = 1 / (CoV^2 * N);
        err_FORM = abs(Pf_FORM - Pf_MC) / Pf_MC;
        err_SORM = abs(Pf_SORM - Pf_MC) / Pf_MC;
        failures = Pf_MC * N;
    end

    %% Store results
    results(i).c          = c;
    results(i).N          = N;
    results(i).Pf_MC      = Pf_MC;
    results(i).Pf_FORM    = Pf_FORM;
    results(i).Pf_SORM    = Pf_SORM;
    results(i).CoV        = CoV;
    results(i).efficiency = efficiency;
    results(i).err_FORM   = err_FORM;
    results(i).err_SORM   = err_SORM;
    results(i).failures   = failures;

    %% Console output
    fprintf('MC:   Pf = %.3e | CoV = %.3f | failures ≈ %.0f\n', Pf_MC, CoV, failures);
    fprintf('FORM: Pf = %.3e | err = %.3f\n', Pf_FORM, err_FORM);
    fprintf('SORM: Pf = %.3e | err = %.3f\n\n', Pf_SORM, err_SORM);

end

end