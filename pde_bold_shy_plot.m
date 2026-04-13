% ==========================================================
% User-defined inputs for this run
% ==========================================================
toxin_level = 90;
weak_strong = 1; % weak 1, strong 2

% Import data
X  = importdata("/Users/tianxuwang/Documents/Research/bold-shy pde/code/X.txt");
Y1 = importdata("/Users/tianxuwang/Documents/Research/bold-shy pde/code/Y1.txt");
Y2 = importdata("/Users/tianxuwang/Documents/Research/bold-shy pde/code/Y2.txt");
Z  = importdata("/Users/tianxuwang/Documents/Research/bold-shy pde/code/Z.txt");

% Set parameters
A = 100;
T = 3000;
Nt = size(X, 1);
Nx = size(X, 2);

xmesh = linspace(-A, A, Nx);
tspan = linspace(0, T, Nt);
dx = xmesh(2) - xmesh(1);

% Last quarter of the time period
idx_t = floor(3*Nt/4) + 1 : Nt;

% 1) Space-time average density over the last quarter
avg_pop1 = mean(Y1(idx_t, :), 'all');
avg_pop2 = mean(Y2(idx_t, :), 'all');
avg_pop_total = mean(Y1(idx_t, :) + Y2(idx_t, :), 'all');

% 2) Sum over space, then average over time, over the last quarter
pop1_space_int = sum(Y1(idx_t, :), 2) * dx;
pop2_space_int = sum(Y2(idx_t, :), 2) * dx;
pop_total_space_int = sum(Y1(idx_t, :) + Y2(idx_t, :), 2) * dx;

sum_avg_pop1 = mean(pop1_space_int);
sum_avg_pop2 = mean(pop2_space_int);
sum_avg_pop_total = mean(pop_total_space_int);

% Display results
disp([toxin_level, weak_strong, avg_pop1, avg_pop2, avg_pop_total, ...
      sum_avg_pop1, sum_avg_pop2, sum_avg_pop_total])

% ==========================================================
% Plot figures
% ==========================================================
tt = tspan;
beg = 1;

% Figure for X
figure('Name', 'X')
surf(xmesh, tt(beg:Nt), X(beg:Nt, :), 'edgecolor', 'none')
ylim([0 T])
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$t$', 'Interpreter', 'latex')
zlabel('$X$', 'Interpreter', 'latex')
view([172 67])
ax = gca;
ax.FontSize = 20;

% Figure for Y1
figure('Name', 'Y1')
surf(xmesh, tt(beg:Nt), Y1(beg:Nt, :), 'edgecolor', 'none')
ylim([0 T])
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$t$', 'Interpreter', 'latex')
zlabel('$Y_1$', 'Interpreter', 'latex')
view([172 67])
ax = gca;
ax.FontSize = 20;

% Figure for Y2
figure('Name', 'Y2')
surf(xmesh, tt(beg:Nt), Y2(beg:Nt, :), 'edgecolor', 'none')
ylim([0 T])
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$t$', 'Interpreter', 'latex')
zlabel('$Y_2$', 'Interpreter', 'latex')
view([172 67])
ax = gca;
ax.FontSize = 20;

% Figure for Z
figure('Name', 'Z')
surf(xmesh, tt(beg:Nt), Z(beg:Nt, :), 'edgecolor', 'none')
ylim([0 T])
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$t$', 'Interpreter', 'latex')
zlabel('$Z$', 'Interpreter', 'latex')
view([172 67])
ax = gca;
ax.FontSize = 20;

% ==========================================================
% Save one row per run into a CSV file
% ==========================================================
outfile = "/Users/tianxuwang/Documents/Research/bold-shy pde/code/population_summary.csv";

newRow = table( ...
    toxin_level, weak_strong, ...
    avg_pop1, avg_pop2, avg_pop_total, ...
    sum_avg_pop1, sum_avg_pop2, sum_avg_pop_total, ...
    'VariableNames', {'toxin_level', 'weak_strong', ...
                      'avg_pop1', 'avg_pop2', 'avg_pop_total', ...
                      'sum_avg_pop1', 'sum_avg_pop2', 'sum_avg_pop_total'} );

if isfile(outfile)
    writetable(newRow, outfile, 'WriteMode', 'append');
else
    writetable(newRow, outfile);
end