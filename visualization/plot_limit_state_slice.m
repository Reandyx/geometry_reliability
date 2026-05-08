function plot_limit_state_slice(problem, u_star, idx1, idx2)
%PLOT_LIMIT_STATE_SLICE 2D slice of high-dimensional limit-state
%
% Inputs:
%   problem  - problem struct
%   u_star   - design point (n-dim vector)
%   idx1     - first dimension to vary
%   idx2     - second dimension to vary

    n = length(u_star);

    x = linspace(-4, 4, 200);
    [X1, X2] = meshgrid(x, x);

    g_vals = zeros(size(X1));

    for i = 1:numel(X1)

        u = u_star; % start at design point

        % vary selected dimensions
        u(idx1) = X1(i);
        u(idx2) = X2(i);

        g_vals(i) = evaluate_limit_state(problem, u);
    end

    G = reshape(g_vals, size(X1));

    % =========================
    % Plot
    % =========================
    figure;
    hold on;

    % Failure region
    contourf(X1, X2, G, [-1e6 0], 'LineStyle', 'none');
    colormap([1 0.8 0.8]);

    % Limit-state boundary
    contour(X1, X2, G, [0 0], 'k', 'LineWidth', 2);

    % Design point projection
    plot(u_star(idx1), u_star(idx2), 'ro', ...
        'MarkerSize', 8, 'LineWidth', 2);

    xlabel(sprintf('u_%d', idx1));
    ylabel(sprintf('u_%d', idx2));
    title('Limit-State Slice');

    axis equal;
    grid on;
end