function method = select_method_from_geometry(class, gamma, beta, topology_flag)

% ==============================
% Priority 1 — TOPOLOGY
% ==============================
if topology_flag == 1
    method = 'MCS';
    return;
end

% ==============================
% Priority 2 — RARE EVENT (FEASIBILITY)
% ==============================
if beta >= 5
    method = 'IS';
    return;
end

% ==============================
% Priority 3 — CURVATURE (ACCURACY)
% ==============================
switch class

    case 'linear'
        method = 'FORM';

    case 'moderate_curvature'
        method = 'SORM';

    case 'high_curvature'
        method = 'SORM';

    case 'rare_event'
        % DO NOT override curvature unless beta triggers IS
        if gamma < 0.1
            method = 'FORM';
        else
            method = 'SORM';
        end

    otherwise
        method = 'SORM';

end

end