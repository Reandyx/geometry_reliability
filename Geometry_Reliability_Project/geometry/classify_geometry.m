function class = classify_geometry(N, kappa, topology_flag, Pf)
%CLASSIFY_GEOMETRY Classifies reliability problem geometry

    % -----------------------------
    % INPUT VALIDATION
    % -----------------------------
    if nargin < 4
        Pf = 1; % default (non-rare)
    end

    if isempty(N) || ~isnumeric(N)
        error('N must be a numeric scalar');
    end

    if isempty(topology_flag)
        topology_flag = 0;
    end

    % -----------------------------
    % CLASSIFICATION (PRIORITY ORDER)
    % -----------------------------

    % TOPOLOGY (strongest effect)
    if topology_flag == 1
        class = 'disconnected';
        return;
    end

    % RARE EVENT (sampling breakdown)
    if Pf < 1e-4
        class = 'rare_event';
        return;
    end

    % GEOMETRY (curvature-driven)
    if N < 0.6
        class = 'linear';
    elseif N < 0.75
        class = 'moderate_curvature';
    else
        class = 'high_curvature';
    end

end