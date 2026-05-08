function plot_mpp_failure_map(problem, multi_result)

% =====================================
% AUTOLOAD DEFAULT CASE
% =====================================
if nargin < 2

    data = load('results/multipoint/raw/data.mat');

    results = data.results;

    % Default showcase case
    c_index = 4;

    c = results(c_index).c;

    b = calibrate_b_disconnected_c(c, 0.15);

    problem = get_problem_disconnected_c(b, c);

    multi_result = results(c_index).MULTI;

end
% =====================================
% Grid
% =====================================
x1 = linspace(-4,4,300);

x2 = linspace(-4,4,300);

[X1,X2] = meshgrid(x1,x2);

G = zeros(size(X1));

% =====================================
% Evaluate limit-state
% =====================================
for i = 1:size(X1,1)

    for j = 1:size(X1,2)

        U = [X1(i,j), X2(i,j)];

        G(i,j) = ...
            evaluate_limit_state_u(problem, U);

    end

end

% =====================================
% Plot failure surface
% =====================================
figure;

contourf(X1, X2, G, 40, ...
    'LineColor','none');

hold on;

% Failure boundary
contour(X1, X2, G, [0 0], ...
    'k', 'LineWidth', 2);

colorbar;

xlabel('u_1');

ylabel('u_2');

title('Detected MPPs over Failure Surface');

grid on;

axis equal;

% =====================================
% Plot MPPs
% =====================================
branches = multi_result.branch_results;

for i = 1:length(branches)

    U_star = branches(i).U_star;

    plot(U_star(1), U_star(2), ...
        'ro', ...
        'MarkerSize', 10, ...
        'LineWidth', 2);

    text(U_star(1)+0.1, ...
         U_star(2), ...
         sprintf('MPP %d', i), ...
         'FontSize', 10, ...
         'Color', 'w');

end

% =====================================
% Origin
% =====================================
plot(0,0,'wx', ...
    'MarkerSize',12, ...
    'LineWidth',2);

legend('Failure region', ...
       'Failure boundary', ...
       'MPP');

end