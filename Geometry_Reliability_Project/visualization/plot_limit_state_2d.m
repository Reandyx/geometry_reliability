function plot_limit_state_2d(problem, u_star)
%PLOT_LIMIT_STATE_2D Visualizes limit-state function in 2D
%
% Inputs:
%   problem  - struct with evaluate_limit_state
%   u_star   - FORM design point [u1, u2]

    % Grid
    x = linspace(-4, 4, 200);
    [X1, X2] = meshgrid(x, x);

    g_vals = zeros(size(X1));

    % Evaluate g over grid
    for i = 1:numel(X1)
        u = [X1(i), X2(i)];
        g_vals(i) = evaluate_limit_state(problem, u);
    end

    G = reshape(g_vals, size(X1));

    % =========================
    % Plot
    % =========================
    figure;
    hold on;

    % --- Failure region (g < 0) ---
    contourf(X1, X2, reshape(g_vals, size(X1)) < 0, ...
    [1 1], 'FaceAlpha', 0.3, 'LineStyle', 'none');
    colormap([1 0.8 0.8]); % light red

    % --- Limit-state boundary g=0 ---
    contour(X1, X2, G, [0 0], 'k', 'LineWidth', 2);

    % --- Design point ---
    plot(u_star(1), u_star(2), 'ro', ...
        'MarkerSize', 8, 'LineWidth', 2);

    % --- Origin ---
    plot(0, 0, 'kx', 'LineWidth', 2);

    xlabel('u_1');
    ylabel('u_2');
    title('Limit-State Surface and Failure Region');

    legend('Failure region', 'g=0', 'Design point', 'Origin');

    axis equal;
    grid on;
end