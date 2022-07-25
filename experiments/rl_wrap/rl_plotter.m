classdef rl_plotter < plotter_hy_interface
    %RL_PLOTTER Plotter for right-left wrap experiment
    %   Detailed explanation goes here
    
    properties        
        R0 = 0.2;
        C0 = [-0.75; 0.75];
        R1 = 1;
        Cu = [0;0];
    end
    
    methods
        function obj = rl_plotter(osm, osd, R0, C0, Cu)
            %CUBE_PLOTTER Construct an instance of this class
            %   Detailed explanation goes here
%             obj.Property1 = inputArg1 + inputArg2;
            obj@plotter_hy_interface(osm, osd)
            
            if nargin > 2
%                 obj.R = R;
                obj.R0 = R0;
                obj.C0 = C0;
                
            end

            if nargin > 4
                obj.Cu = Cu;
            end
        end
        
        function [F, ax] = rl_plot(obj,p_est)
            %plot trajectories of the system in the unit box [-R,R]
            
            F = figure(800);
            clf
            hold on
            

            loc_names = {''};
            for i = 1:length(obj.osd.locations) %i=1:2
        %         subplot(1, 2, i)
        %         hold on
                xlabel('x_1')
                ylabel('x_2')
                
                title( loc_names, 'FontSize', 14)


%                 [X, Y, Z] = sphere(30);
%                 surf( obj.R0*X, obj.R0*Y 'FaceColor', 0.5*[1,1,1], 'FaceAlpha', obj.falpha, 'edgecolor', 'none');

        %         title(['Location ', num2str(i)], 'FontSize', 14)

                %trajectory
                %exploit the symmetry structure of dynamics to plot twice
                %as many trajectories
                for j = 1:length(obj.osd.locations{i})
                    traj_curr = obj.osd.locations{i}{j};
                    plot( traj_curr.x(:, 1), traj_curr.x(:, 2), 'c')

                end 

            end    
            xlim([-obj.R1, obj.R1])
            ylim([-obj.R1, obj.R1])
            axis square

            %plot circles
            theta = linspace(0, 2*pi, 300);
            circ = [cos(theta); sin(theta)];
            X0 = obj.R0*circ + obj.C0;
            plot(X0(1, :), X0(2, :), 'k', 'linewidth', 3)

            
            
            %TODO: add visualization of p_est patches
            
            
            
            if nargin == 2
                Rp = sqrt(-p_est);
                Xp = Rp * circ + obj.Cu;
                plot(Xp(1, :), Xp(2, :), ':r', 'linewidth', 3)
            end
%             linkprop([ax2; ax1], {'XLim', 'YLim', 'ZLim', 'CameraPosition', 'CameraUpVector', 'CameraViewAngle', 'CameraTarget'});
        end
    end
end

