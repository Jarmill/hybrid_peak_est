classdef cube_plotter < plotter_hy_interface
    %CUBE_PLOTTER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        R0 = 0.2;
        R1 = 1;
        falpha = 0.4;
    end
    
    methods
        function obj = cube_plotter(osm, osd, R0, R1)
            %CUBE_PLOTTER Construct an instance of this class
            %   Detailed explanation goes here
%             obj.Property1 = inputArg1 + inputArg2;
            obj@plotter_hy_interface(osm, osd)
            
            if nargin > 2
                obj.R0 = R0;
                obj.R1 = R1;
            end
        end
        
        function [F, ax1, ax2] = cube_plot(obj,p_est)
            %plot trajectories of the system in the unit cube [-1,1]^3
            
            F = figure(1);
            clf

            loc_names = {'No Control', 'Controlled'};
%             ax1 = subplot(1,2,1);
            ax1 = subplot(2,1,1);
            hold on
            axis off
pbaspect([1,1,1])
%             ax2 = subplot(1,2,2);
            ax2 = subplot(2,1,2);
            hold on
            axlist = [ax1, ax2];
            axis off
pbaspect([1,1,1])

            for i = 1:length(obj.osd.locations) %i=1:2
        %         subplot(1, 2, i)
        %         hold on
                xlabel(axlist(i), 'x_1')
                ylabel(axlist(i), 'x_2')
                zlabel(axlist(i), 'x_3')
                

%                 title(axlist(i), loc_names{i}, 'FontSize', 16)
                [X, Y, Z] = sphere(30);
                surf(axlist(i), obj.R0*X, obj.R0*Y, obj.R0*Z, 'FaceColor', 0.5*[1,1,1], 'FaceAlpha', obj.falpha, 'edgecolor', 'none');

        %         title(['Location ', num2str(i)], 'FontSize', 14)

                %trajectory
                %exploit the symmetry structure of dynamics to plot twice
                %as many trajectories
                for j = 1:length(obj.osd.locations{i})
                    traj_curr = obj.osd.locations{i}{j};
                    plot3(axlist(i), traj_curr.x(:, 1), traj_curr.x(:, 2), traj_curr.x(:, 3), 'c')
                    plot3(axlist(i),-traj_curr.x(:, 1),-traj_curr.x(:, 2),-traj_curr.x(:, 3), 'c')
                    
                    
%                     if j < length(obj.osd.locations{i})
%                         scatter3(axlist(i), traj_curr.x(end, 1), traj_curr.x(end, 2), traj_curr.x(end, 3), 100, 'ok')
%                         scatter3(axlist(i),-traj_curr.x(end, 1),-traj_curr.x(end, 2),-traj_curr.x(end, 3), 100, 'ok')
%                     end
%                     if (i==1) && (j==1)
%                         scatter3(axlist(i), traj_curr.x(1, 1), traj_curr.x(1, 2), traj_curr.x(1, 3), 400, 'ok')
%                         scatter3(axlist(i),-traj_curr.x(1, 1),-traj_curr.x(1, 2),-traj_curr.x(1, 3), 400, 'ok')
%                     else
%                         scatter3(axlist(i), traj_curr.x(1, 1), traj_curr.x(1, 2), traj_curr.x(1, 3), 100, 'ok')
%                         scatter3(axlist(i),-traj_curr.x(1, 1),-traj_curr.x(1, 2),-traj_curr.x(1, 3), 100, 'ok')
%                     end



                end 
        %         axis square
                hlink = linkprop([ax1,ax2],{'CameraPosition','CameraUpVector', 'CameraTarget'});
                setappdata(gcf, 'theLink', hlink)
                view(3)

%                 xlim(1.1*p_est*[-1;1])


            end    
            
            %TODO: add visualization of p_est patches
            
            
            
            if nargin == 2
                yl = ylim;
                zl = zlim;
                for i = 1:2
                    xp = ones(1, 5);
                    yp = yl([1 2 2 1 1]);
                    zp = zl([1 1 2 2 1]);
                    if i==2
                    patch(axlist(i),p_est*xp, yp, zp, 'r','FaceAlpha', obj.falpha, 'edgecolor', 'none')
                    patch(axlist(i),-p_est*xp, yp, zp, 'r','FaceAlpha', obj.falpha, 'edgecolor', 'none')
                    end
                    
                    
                end
            end
            linkprop([ax2; ax1], {'XLim', 'YLim', 'ZLim', 'CameraPosition', 'CameraUpVector', 'CameraViewAngle', 'CameraTarget'});
        end
    end
end

