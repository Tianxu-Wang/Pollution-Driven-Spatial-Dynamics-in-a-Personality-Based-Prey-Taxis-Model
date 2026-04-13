X  = importdata("/Users/tianxuwang/Documents/Research/bold-shy pde/code/X.txt");
Y1 = importdata("/Users/tianxuwang/Documents/Research/bold-shy pde/code/Y1.txt");
Y2 = importdata("/Users/tianxuwang/Documents/Research/bold-shy pde/code/Y2.txt");
Z  = importdata("/Users/tianxuwang/Documents/Research/bold-shy pde/code/Z.txt");

Nt = size(X, 1);
Nx = 201;
Ny = 201;
X  = reshape(X,  [Nt Nx Ny]);
Y1 = reshape(Y1, [Nt Nx Ny]);
Y2 = reshape(Y2, [Nt Nx Ny]);
Z  = reshape(Z,  [Nt Nx Ny]);
A = 100;
T = 500;
xmesh = linspace(-A, A, Nx);
ymesh = linspace(-A, A, Ny);
tspan = linspace(0, T, Nt);




% numFrames = 20;
% idx = round(linspace(1, Nt, numFrames));


% Create folders for .fig and .png files
figFolder = 'snapshots_fig_T30_threecenter';
pngFolder = 'snapshots_png_T30_threecenter';
if ~exist(figFolder, 'dir'), mkdir(figFolder); end
if ~exist(pngFolder, 'dir'), mkdir(pngFolder); end

% Choose 20 equally spaced time indices
numFrames = 100;
idx = round(linspace(1, Nt, numFrames));

% Loop over selected time indices
% for k = 1:length(idx)
% % it = idx(k);
for k = 500:501
    it = k;
    t_val = tspan(it);

    % ----- X -----
    fig1 = figure('Visible', 'off');
    surf(xmesh, ymesh, reshape(X(it, :, :), Nx, Ny), 'edgecolor', 'none')
    xlabel('$x$','Interpreter','latex')
    ylabel('$y$','Interpreter','latex')
    zlabel('X(t,x,y)')
    ax = gca; ax.FontSize = 20;
    % clim([0 1.2]);
    colormap("turbo");
    colorbar;
    view(2)
    % savefig(fig1, fullfile(figFolder, sprintf('X_t%.2f_T10.fig', t_val)));
    exportgraphics(fig1, fullfile(pngFolder, sprintf('X_t%.2f_T10.png', t_val)), 'Resolution', 300);
     close(fig1);

    % ----- Y1 -----
    fig2 = figure('Visible', 'off');
    surf(xmesh, ymesh, reshape(Y1(it, :, :), Nx, Ny), 'edgecolor', 'none')
    xlabel('$x$','Interpreter','latex')
    ylabel('$y$','Interpreter','latex')
    zlabel('Y1(t,x,y)')
    ax = gca; ax.FontSize = 20;
    colormap("turbo");
    colorbar;
    view(2)
      % clim([0 1]);
    % savefig(fig2, fullfile(figFolder, sprintf('Y1_t%.2f_T10.fig', t_val)));
    exportgraphics(fig2, fullfile(pngFolder, sprintf('Y1_t%.2f_T0.png', t_val)), 'Resolution', 300);
     close(fig2);

    % % ----- Y2 -----
    fig3 = figure('Visible', 'off');
    surf(xmesh, ymesh, reshape(Y2(it, :, :), Nx, Ny), 'edgecolor', 'none')
    xlabel('$x$','Interpreter','latex')
    ylabel('$y$','Interpreter','latex')
    zlabel('Y2(t,x,y)')
    ax = gca; ax.FontSize = 20;
    colormap("turbo");
    colorbar;
      % clim([0 0.1]);
    view(2)
    % savefig(fig3, fullfile(figFolder, sprintf('Y2_t%.2f_T10.fig', t_val)));
    exportgraphics(fig3, fullfile(pngFolder, sprintf('Y2_t%.2f_T0.png', t_val)), 'Resolution', 300);
     close(fig3);

    % % ----- Z -----
    fig4 = figure('Visible', 'off');
    surf(xmesh, ymesh, reshape(Z(it, :, :), Nx, Ny), 'edgecolor', 'none')
    xlabel('$x$','Interpreter','latex')
    ylabel('$y$','Interpreter','latex')
    zlabel('Z(t,x,y)')
    ax = gca; ax.FontSize = 20;
    colormap("turbo");
    colorbar;
    view(2)
    % clim([0 1]);
    % savefig(fig4, fullfile(figFolder, sprintf('Z_t%.2f_T10.fig', t_val)));
    % FIX: changed .fig to .png for the PNG export
    exportgraphics(fig4, fullfile(pngFolder, sprintf('Z_t%.2f_T10.png', t_val)), 'Resolution', 300);
     close(fig4);

    fprintf('Saved snapshot %d of %d (t = %.2f)\n', k, numFrames, t_val);
end
% % Loop only over those indices
% 
% 
% 
% for it = 500:501
% 
%     figure(1)
%     surf(xmesh, ymesh, reshape(X(it, :, :), Nx, Ny), 'edgecolor', 'none')
%     xlabel('$x$','Interpreter','latex')
%     ylabel('$y$','Interpreter','latex')
%     zlabel('X(t,x,y)')
%     title(sprintf('$X$ at $t = %.f$', tspan(it)), 'FontSize', 20,Interpreter='latex')
%     ax = gca;
%     ax.FontSize = 20;
%     % caxis([0 1.2]);
%     colormap("turbo");
%     colorbar;
%     view(2)
% 
%     figure(2)
%     surf(xmesh, ymesh, reshape(Y1(it, :, :), Nx, Ny), 'edgecolor', 'none')
%     xlabel('$x$','Interpreter','latex')
%     ylabel('$y$','Interpreter','latex')
%     zlabel('Y1(t,x,y)')
%     title(sprintf('$Y_1$ at $t = %.f$', tspan(it)), 'FontSize', 20, Interpreter='latex')
%     ax = gca;
%     ax.FontSize = 20;
%     colormap("turbo");
%     colorbar;
%     % clim([0 2]);
%     view(2)
% 
%     figure(3)
%     surf(xmesh, ymesh, reshape(Y2(it, :, :), Nx, Ny), 'edgecolor', 'none')
%     xlabel('$x$','Interpreter','latex')
%     ylabel('$y$','Interpreter','latex')
%     zlabel('Y2(t,x,y)')
%     title(sprintf('$Y_2$ at $t = %.0f$', tspan(it)), 'FontSize', 20, Interpreter='latex')
%     ax = gca;
%     ax.FontSize = 20;
%     colormap("turbo");
%     colorbar;
%     % clim([0 0.5]);
%     view(2)
% 
%     figure(4)
%     surf(xmesh, ymesh, reshape(Z(it, :, :), Nx, Ny), 'edgecolor', 'none')
%     xlabel('$x$','Interpreter','latex')
%     ylabel('$y$','Interpreter','latex')
%     zlabel('Z(t,x,y)')
%     title(sprintf('$Z$ at $t = %.0f$', tspan(it)), 'FontSize', 20, Interpreter='latex')
%     ax = gca;
%     ax.FontSize = 20;
%     colormap("turbo");
%     colorbar;
%     % clim([0 1]);
%     view(2)
% 
%     pause(0.1);
% end
