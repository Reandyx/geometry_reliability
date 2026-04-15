function results = run_curvature_study(model_type)

    if nargin < 1
        model_type = 'global';
    end

    a_values = [0.5, 1, 2, 5, 10];
    Pf_target = 0.2;
    N_mc = 1e6;

    n_cases = length(a_values);

    results(n_cases) = struct( ...
        'a', [], ...
        'model', [], ...
        'Pf_MC', [], ...
        'Pf_FORM', [], ...
        'Pf_SORM', [], ...
        'error_FORM', [], ...
        'error_SORM', [], ...
        'improvement', [], ...
        'beta_FORM', [], ...
        'sigma_MC', [], ...
        'bias', [], ...
        'noise_ratio', [], ...
        'U_star', [], ...
        'grad_g', [], ...
        'H', [], ...
        'kappa_global', [], ...
        'kappa_principal', [], ...
        'N', []);

    fprintf('\n============================\n');
    fprintf('   GEOMETRY + SORM STUDY (M6)\n');
    fprintf('============================\n');

    for i = 1:n_cases

        a = a_values(i);

        fprintf('\n--- CASE a = %.2f ---\n', a);

        % =========================
        % Problem setup
        % =========================
        b = calibrate_b(a, Pf_target, model_type);
        fprintf('Calibrated b = %.4f\n', b);

        problem = build_problem(a, b, model_type);

        % =========================
        % Monte Carlo
        % =========================
        reset_eval_count();
        res_mc = run_mcs(problem, N_mc);

        Pf_mc = res_mc.Pf;
        sigma_mc = sqrt(Pf_mc * (1 - Pf_mc) / N_mc);

        fprintf('MCS  Pf = %.6e ± %.2e\n', Pf_mc, sigma_mc);

        % =========================
        % FORM
        % =========================
        reset_eval_count();
        res_form = run_form(problem);

        Pf_form = res_form.Pf;
        bias = Pf_form - Pf_mc;

        fprintf('FORM Pf = %.6e\n', Pf_form);

        % =========================
        % SORM
        % =========================
        res_sorm = run_sorm(problem);
        Pf_sorm = res_sorm.Pf;

        fprintf('SORM Pf = %.6e\n', Pf_sorm);

        % =========================
        % Geometry (M5 + M6)
        % =========================
        U_star = res_form.U_star;
        grad_g = res_form.grad_g;

        H = compute_hessian_u(problem, U_star);

        kappa_principal = compute_principal_curvatures(H, grad_g);
        
        if isempty(kappa_principal)
            kappa_global = 0;
        else
            kappa_global = norm(kappa_principal);
        end
        
        N = kappa_global;

        fprintf('kappa_global = %.4e\n', kappa_global);
        fprintf('N index      = %.4e\n', N);
        
        fprintf('Principal curvatures: ');
        fprintf('%.3e ', kappa_principal);
        fprintf('\n');

        % =========================
        % Errors
        % =========================
        if Pf_mc == 0
            error_FORM = NaN;
            error_SORM = NaN;
        else
            error_FORM = abs(Pf_form - Pf_mc) / Pf_mc;
            error_SORM = abs(Pf_sorm - Pf_mc) / Pf_mc;
        end

        improvement = error_FORM - error_SORM;

        % MC noise awareness
        if sigma_mc > 0
            noise_ratio = abs(Pf_form - Pf_mc) / sigma_mc;
        else
            noise_ratio = NaN;
        end

        fprintf('Error FORM = %.6f\n', error_FORM);
        fprintf('Error SORM = %.6f\n', error_SORM);
        fprintf('Improvement = %.6f\n', improvement);
        fprintf('NoiseR  = %.4f\n', noise_ratio);
        fprintf('Beta    = %.4f\n', res_form.beta);

        % =========================
        % Store
        % =========================
        results(i).a = a;
        results(i).model = model_type;

        results(i).Pf_MC = Pf_mc;
        results(i).Pf_FORM = Pf_form;
        results(i).Pf_SORM = Pf_sorm;

        results(i).error_FORM = error_FORM;
        results(i).error_SORM = error_SORM;
        results(i).improvement = improvement;

        results(i).beta_FORM = res_form.beta;

        results(i).sigma_MC = sigma_mc;
        results(i).bias = bias;
        results(i).noise_ratio = noise_ratio;

        results(i).U_star = U_star;
        results(i).grad_g = grad_g;

        results(i).H = H;

        results(i).kappa_global = kappa_global;
        results(i).kappa_principal = kappa_principal;
        results(i).N = N;
        
        results(i).problem = problem;
        
    end

    fprintf('\n============================\n');
    fprintf('   STUDY COMPLETE\n');
    fprintf('============================\n');

end