%% Import data
X  = importdata("/Users/tianxuwang/Documents/Research/bold-shy pde/code/X.txt");
Y1 = importdata("/Users/tianxuwang/Documents/Research/bold-shy pde/code/Y1.txt");
Y2 = importdata("/Users/tianxuwang/Documents/Research/bold-shy pde/code/Y2.txt");
Z  = importdata("/Users/tianxuwang/Documents/Research/bold-shy pde/code/Z.txt");

%% Parameters

A = 100;
T = 3000;
Tview = 3000;

Nt = size(X, 1);
Nx = size(X, 2);

xmesh = linspace(-A, A, Nx);
tspan = linspace(0, T, Nt);

%% Colors from your figure
cX  = [0.0000 0.4470 0.7410];
cY1 = [0.9294 0.6941 0.1255];
cY2 = [0.8510 0.3255 0.0980];
cZ  = [0.4667 0.6745 0.1882];

%% Restrict data to 0 <= t <= 500
idx500 = find(tspan <= Tview);
tspan500 = tspan(idx500);

X500  = X(idx500, :);
Y1500 = Y1(idx500, :);
Y2500 = Y2(idx500, :);
Z500  = Z(idx500, :);

%% Averages over the last 201 time steps within 0-500
tailLen = min(201, length(idx500));
tailIdxLocal = length(idx500) - tailLen + 1 : length(idx500);

avg_pop1      = mean(Y1500(tailIdxLocal, :), 'all');
avg_pop2      = mean(Y2500(tailIdxLocal, :), 'all');
avg_pop_total = mean(Y1500(tailIdxLocal, :) + Y2500(tailIdxLocal, :), 'all');

fprintf('avg_pop1      = %.6f\n', avg_pop1);
fprintf('avg_pop2      = %.6f\n', avg_pop2);
fprintf('avg_pop_total = %.6f\n', avg_pop_total);

%% Profile at t = 500
Xeq  = X500(end, :);
Y1eq = Y1500(end, :);
Y2eq = Y2500(end, :);
Zeq  = Z500(end, :);

figure('Name', 'Profiles at t = 500', 'Color', 'w');
plot(xmesh, Xeq,  'LineWidth', 1.5, 'Color', cX); hold on;
plot(xmesh, Y1eq, 'LineWidth', 1.5, 'Color', cY1);
plot(xmesh, Y2eq, 'LineWidth', 1.5, 'Color', cY2);
plot(xmesh, Zeq,  'LineWidth', 1.5, 'Color', cZ);
set(gca, 'FontSize', 24);
hold off;

xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 24);
ylabel('Density', 'Interpreter', 'latex', 'FontSize', 24);
legend({'$X$', '$Y_1$', '$Y_2$', '$Z$'}, ...
    'Interpreter', 'latex', 'FontSize', 14, 'Location', 'best');
title(sprintf('Spatial profiles at t \\approx %.2f', tspan500(end)), 'FontSize', 24);
set(gca, 'FontSize', 16, 'LineWidth', 1.2);
xlim([xmesh(1), xmesh(end)]);
grid on;
box on;

%% Time series at the center point x = 0
[~, centerIndex] = min(abs(xmesh));

figure('Name', 'Center point dynamics (0 to 500)', 'Color', 'w');
plot(tspan500, X500(:, centerIndex),  'LineWidth', 1.5, 'Color', cX); hold on;
plot(tspan500, Y1500(:, centerIndex), 'LineWidth', 1.5, 'Color', cY1);
plot(tspan500, Y2500(:, centerIndex), 'LineWidth', 1.5, 'Color', cY2);
plot(tspan500, Z500(:, centerIndex),  'LineWidth', 1.5, 'Color', cZ);
set(gca, 'FontSize', 24);
hold off;

xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 24);
ylabel('Density at $x=0$', 'Interpreter', 'latex', 'FontSize', 24);
legend({'$X(0,t)$', '$Y_1(0,t)$', '$Y_2(0,t)$', '$Z(0,t)$'}, ...
    'Interpreter', 'latex', 'FontSize', 14, 'Location', 'best');
title('Dynamics at the center point for 0 \leq t \leq 500', 'FontSize', 24);
set(gca, 'FontSize', 16, 'LineWidth', 1.2);
xlim([0, Tview]);
grid on;
box on;

% Surface plots only for 0 to 500
figure('Name', 'Surface plots (0 to 500)', 'Color', 'w');
tiledlayout(2,2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
surf(xmesh, tspan500, X500, 'EdgeColor', 'none');
colormap(gca, parula);
view(172, 67);
% view(45, 30)
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$t$', 'Interpreter', 'latex');
zlabel('$X$', 'Interpreter', 'latex');
title('$X$', 'Interpreter', 'latex');
set(gca, 'FontSize', 14);
ylim([0, Tview]);
colorbar;

nexttile;
surf(xmesh, tspan500, Y1500, 'EdgeColor', 'none');
colormap(gca, parula);
view(172, 67);
% view(45, 30)
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$t$', 'Interpreter', 'latex');
zlabel('$Y_1$', 'Interpreter', 'latex');
title('$Y_1$', 'Interpreter', 'latex');
set(gca, 'FontSize', 14);
ylim([0, Tview]);
colorbar;

nexttile;
surf(xmesh, tspan500, Y2500, 'EdgeColor', 'none');
colormap(gca, parula);
view(172, 67);
% view(45, 30)
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$t$', 'Interpreter', 'latex');
zlabel('$Y_2$', 'Interpreter', 'latex');
title('$Y_2$', 'Interpreter', 'latex');
set(gca, 'FontSize', 14);
ylim([0, Tview]);
colorbar;

nexttile;
surf(xmesh, tspan500, Z500, 'EdgeColor', 'none');
colormap(gca, parula);
 view(172, 67);
% view(45, 30)
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$t$', 'Interpreter', 'latex');
zlabel('$Z$', 'Interpreter', 'latex');
title('$Z$', 'Interpreter', 'latex');
set(gca, 'FontSize', 14);
ylim([0, Tview]);
colorbar;