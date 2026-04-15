function plot_disconnected_domain(problem, U_star)

    % grid in X-space
    x1 = linspace(-5,5,300);
    x2 = linspace(-5,5,300);

    [X1, X2] = meshgrid(x1, x2);

    X = [X1(:), X2(:)];

    g = problem.gfun(X);
    G = reshape(g, size(X1));

    figure;
    contour(X1, X2, G, [0 0], 'LineWidth', 2);
    hold on;

    % failure region shading
    contourf(X1, X2, G, [-10 0]);
    colormap([1 0.8 0.8]);

    % design point (convert U -> X)
    X_star = u_to_x(U_star', problem);

    plot(X_star(1), X_star(2), 'ro', 'MarkerSize', 10, 'LineWidth', 2);

    xlabel('x1');
    ylabel('x2');
    title('Disconnected Failure Domains + FORM Design Point');
    grid on;

end