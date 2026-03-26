function result = create_result_struct()

    % Core fields (REQUIRED everywhere)
    result.method     = '';
    result.Pf         = NaN;
    result.beta       = NaN;
    result.neval      = 0;
    result.runtime    = 0;
    result.converged  = false;

    % FORM / SORM shared
    result.U_star     = [];
    result.alpha      = [];
    result.grad_g     = [];
    result.kappa      = [];

    % Debug / optional
    result.history    = [];
    result.g_eval     = [];

end