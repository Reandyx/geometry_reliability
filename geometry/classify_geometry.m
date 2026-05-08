function class = classify_geometry(gamma, topology_flag, Pf)
% CLASSIFY_GEOMETRY
% Geometry classification based on canonical gamma metric
%
% gamma = beta * max |kappa_i|

    if nargin < 3
        Pf = 1;
    end

    if isempty(topology_flag)
        topology_flag = 0;
    end

    % ==============================
    % PRIORITY 1 — TOPOLOGY
    % ==============================
    if topology_flag == 1
        class = 'disconnected';
        return;
    end

    % ==============================
    % PRIORITY 2 — RARE EVENT
    % ==============================
    if Pf < 1e-6
        class = 'rare_event';
        return;
    end

    % ==============================
    % PRIORITY 3 — CURVATURE (GAMMA)
    % ==============================
    if gamma < 0.1
        class = 'linear';
    elseif gamma < 1.0
        class = 'moderate_curvature';
    else
        class = 'high_curvature';
    end

end