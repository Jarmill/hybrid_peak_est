load('traj_stoch_cube_half', 'osd')

param = struct('R0', 0.2, 'R1', 1.5, 'K', 0.25, 'sigma', 0.1, 'dt', 1e-3, 'L', 1.5);

%% plot
%cube_plotter is not playing nicely. I will try and visualize this myself.

            F = figure(100);
            clf

loc_names = {'No Control', 'Controlled'};
tl = tiledlayout(1, 2);
% ax1 = subplot(1,2,1);
ax1 = nexttile;
hold on
axis off
pbaspect([1,1,1])
%     view(3)

% ax2 = subplot(1,2,2);
ax2 = nexttile;
hold on

axlist = [ax1, ax2];
axis off
pbaspect([1,1,1])
%     view(3)
    
    for i = 1:length(osd.locations) %i=1:2
        xlabel(axlist(i), 'x_1')
        ylabel(axlist(i), 'x_2')
        zlabel(axlist(i), 'x_3')

        [X, Y, Z] = sphere(30);
        surf(axlist(i), param.R0*X, param.R0*Y, param.R0*Z, 'FaceColor', 0.5*[1,1,1], 'FaceAlpha', 0.5, 'edgecolor', 'none');

        %trajectory
        %exploit the symmetry structure of dynamics to plot twice
        %as many trajectories
        for j = 1:length(osd.locations{i})
%         for j = 1:2
            traj_curr = osd.locations{i}{j};
            plot3(axlist(i), traj_curr.x(:, 1), traj_curr.x(:, 2), traj_curr.x(:, 3), 'c')
            plot3(axlist(i),-traj_curr.x(:, 1),-traj_curr.x(:, 2),-traj_curr.x(:, 3), 'c')                                        
        end 
        
        linkprop([ax2; ax1], {'XLim', 'YLim', 'ZLim', 'CameraPosition', 'CameraUpVector', 'CameraViewAngle', 'CameraTarget'});
        
        xlim([-1, 1]*param.L)
        ylim([-1, 1]*param.L)
        zlim([-1, 1]*param.L)
    end    
    
        view(3)
   