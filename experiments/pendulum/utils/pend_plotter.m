classdef pend_plotter
    %PEND_PLOTTER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        FS_axis = 14;       %font size of axis
        FS_title = 16;      %font size of title
        G_size = 40;        %size of guard change linespec
%         w_max_est = 0;      %estimated maximum angular velocity
    end
    
    properties(Access=private)
        osm;
    end
    
    methods
        function obj = pend_plotter(osm)
            %PEND_PLOTTER Plot the output of the pendulum
            %   OSM: out_sim_multi
            %   output of pend_manager.pend_sample_multi
%             obj.osm = osm;
%             obj.set_osm(osm);
            if iscell(osm)
                obj.osm = osm;
            else
                obj.osm = {osm};
            end
%             obj.w_max_est = w_max_est;
        end
        
        function obj = set_osm(obj, osm_new)
            obj.osm = osm_new;
        end
        
        function osm_out = get_osm(obj)
            osm_out = obj.osm;
        end
        
        function F = nonneg_jump(obj)
            %NONNEG_JUMP plot nonnegative transitions from the jump            
            F = figure(60);
            clf
            hold on
            title('Jump Nonnegativity', 'FontSize', obj.FS_title)
            ylabel("$v_i(x) - v(R_{i \rightarrow i'} x)$", 'interpreter', 'latex', 'FontSize', obj.FS_axis);
            xlabel('time', 'FontSize', obj.FS_axis)
            for i = 1:length(obj.osm)
                osn = obj.osm{i};
                for j = 1:length(osn.jump)
                    j_curr = osn.jump{j};            
                    stem(j_curr.t, j_curr.nonneg, 'c') 
                end
            end
        end
        
        function F = nonneg_loc(obj)
            %NONNEG_LOC plot nonnegative functions in the locations
            F = figure(61);
            clf
            nonneg_title = {'Initial Value', 'Decrease in Value', 'Cost Proxy'};
            ax_loc = {'$\gamma - v(x)$', '$-L_{f} v(x)$', '$v(x) - p(x)$'};
            for k = 1:3
                subplot(3, 1, k)
                hold on
                xlabel('time', 'FontSize', obj.FS_axis)
                ylabel(ax_loc{k}, 'interpreter', 'latex', 'FontSize', obj.FS_axis);
                title(nonneg_title{k}, 'FontSize', obj.FS_title);   
            end
            for i = 1:length(obj.osm)
                    osn = obj.osm{i};
                for j = 1:length(osn.sim)
                    traj_curr = osn.sim{j};
                    for k = 1:3            
                        subplot(3, 1, k)
                        plot(traj_curr.t, traj_curr.nonneg(:, k), 'c')
                        scatter(traj_curr.t(1), traj_curr.nonneg(1, k), obj.G_size, 'ok')
                        scatter(traj_curr.t(end), traj_curr.nonneg(end, k), obj.G_size, 'ok')
                    end
                end  
            end
        end
        
        function F = time_plot(obj)
            %TIME_PLOT plot states (cos/sin) as a function of time
            F = figure(62);
            clf
            hold on
            for i = 1:length(obj.osm)
                osn = obj.osm{i};
                for j = 1:length(osn.sim)
                    traj_curr = osn.sim{j};
                    plot3(traj_curr.t, traj_curr.xl(:, 2), -traj_curr.xl(:, 1), 'c')    
                    scatter3(traj_curr.t(end), traj_curr.xl(end, 2), -traj_curr.xl(end, 1), obj.G_size, 'ok')
                end
            end
            xlabel('time', 'FontSize', obj.FS_axis);
            ylabel('$sin(\theta)$', 'interpreter', 'latex', 'FontSize', obj.FS_axis);
            zlabel('$-cos(\theta)$', 'interpreter', 'latex', 'FontSize', obj.FS_axis);
            title('Pendulum Swing-up Simulation', 'FontSize',  obj.FS_title);
            view(3)
        end
        
        function F = track_plot(obj)
            %TRACK_PLOT plot quantities (angle, angular velocity, energy)
            F = figure(63);
            clf
            for i = 1:length(obj.osm)
                osn = obj.osm{i};
                for j = 1:length(osn.sim)
                    traj_curr = osn.sim{j};
                    %theta
                    subplot(4, 1, 1)
                    hold on
                    theta_wrap = wrapTo2Pi(traj_curr.x(:, 1));
                    plot(traj_curr.t, theta_wrap, 'c')
%                     scatter(traj_curr.t(1), theta_wrap(1), obj.G_size, 'ok')
                    scatter(traj_curr.t(end), theta_wrap(end), obj.G_size, 'ok')


                    %w
                    subplot(4, 1, 2)
                    hold on
                    plot(traj_curr.t, traj_curr.x(:, 2), 'c')
%                     scatter(traj_curr.t(1), traj_curr.x(1, 2), obj.G_size, 'ok')
                    scatter(traj_curr.t(end), traj_curr.x(end, 2), obj.G_size, 'ok')


                    %energy
                    subplot(4, 1, 3)
                    hold on
                    plot(traj_curr.t, traj_curr.E, 'c')
%                     scatter(traj_curr.t(1), traj_curr.E(1), 40, 'ok')
                    scatter(traj_curr.t(end), traj_curr.E(end), obj.G_size, 'ok')

                    
                                        %energy
                    subplot(4, 1, 4)
                    hold on
                    plot(traj_curr.t, traj_curr.u, 'c')
                    scatter(traj_curr.t(1), traj_curr.u(1), obj.G_size, 'ok')
                    scatter(traj_curr.t(end), traj_curr.u(end), obj.G_size, 'ok')

                    
                end
            end
            %cleanup
            subplot(4, 1, 1)
            ylim([0, 2*pi])
            plot(xlim, [pi, pi], ':k')
            title('Angle', 'FontSize', obj.FS_title)


            subplot(4, 1, 2)
            plot(xlim, [0, 0], ':k')
            title('Angular Velocity', 'FontSize', obj.FS_title)                    
            xlabel('time')

            subplot(4, 1, 3)
            plot(xlim, [1, 1], ':k')
            xlabel('time')
            title('Energy', 'FontSize', obj.FS_title)      
            
            subplot(4, 1, 4)
            plot(xlim, [0, 0], ':k')
            xlabel('time')
            title('Control Effort', 'FontSize', obj.FS_title)  
        end
        
        function F = track_plot_trig(obj)
                        %TRACK_PLOT plot quantities (angle, angular velocity, energy)
            F = figure(64);
            clf
            for i = 1:length(obj.osm)
                osn = obj.osm{i};
                for j = 1:length(osn.sim)
                    traj_curr = osn.sim{j};
                    %cos theta
                    subplot(5, 1, 1)
                    hold on
                    plot(traj_curr.t, traj_curr.xl(:, 1), 'c')
%                     scatter(traj_curr.t(end), traj_curr.xl(end, 2), obj.G_size, 'ok')

                    %sin theta
                    subplot(5, 1, 2)
                    hold on
                    plot(traj_curr.t, traj_curr.xl(:, 2), 'c')
%                     scatter(traj_curr.t(end), traj_curr.xl(end, 2), obj.G_size, 'ok')


                    %w
                    subplot(5, 1, 3)
                    hold on
                    plot(traj_curr.t, traj_curr.x(:, 2), 'c')
%                     scatter(traj_curr.t(1), traj_curr.x(1, 2), obj.G_size, 'ok')
%                     scatter(traj_curr.t(end), traj_curr.x(end, 2), obj.G_size, 'ok')


                    %energy
                    subplot(5, 1, 4)
                    hold on
                    %added 1 to ensure that E > 0
                    plot(traj_curr.t, traj_curr.E+1, 'c')
%                     scatter(traj_curr.t(1), traj_curr.E(1), 40, 'ok')
%                     scatter(traj_curr.t(end), traj_curr.E(end), obj.G_size, 'ok')

                    
                                        %energy
                    subplot(5, 1, 5)
                    hold on
                    plot(traj_curr.t, traj_curr.u, 'c')
%                     scatter(traj_curr.t(1), traj_curr.u(1), obj.G_size, 'ok')
%                     scatter(traj_curr.t(end), traj_curr.u(end), obj.G_size, 'ok')

                    
                end
            end
            %cleanup
            subplot(5, 1, 1)
            ylim([-1, 1])
            plot(xlim, [-1, -1], ':k')            
            title('Cosine of Angle', 'FontSize', obj.FS_title)
            xlabel('time')
            ylabel('$cos(\theta)$', 'interpreter', 'latex')
           
            
            subplot(5, 1, 2)
            ylim([-1, 1])
            plot(xlim, [0, 0], ':k')
            title('Sine of Angle', 'FontSize', obj.FS_title)
            xlabel('time')
            ylabel('$sin(\theta)$', 'interpreter', 'latex')
 
            subplot(5, 1, 3)
            plot(xlim, [0, 0], ':k')
            title('Angular Velocity', 'FontSize', obj.FS_title)                    
            xlabel('time')
            ylabel('$\omega$', 'interpreter', 'latex')
            
            
            subplot(5, 1, 4)
            plot(xlim, [1, 1]+1, ':k')
            xlabel('time')
            title('Energy', 'FontSize', obj.FS_title)      
            ylabel('$E$', 'interpreter', 'latex')
            
            subplot(5, 1, 5)
            plot(xlim, [0, 0], ':k')
            xlabel('time')
            title('Control Effort', 'FontSize', obj.FS_title)  
            ylabel('$u$', 'interpreter', 'latex')
        end
        
        
        function F = state_plot_trig(obj, w_range_0, w_max_est)
            %STATE_PLOT_TRIG plot trajectories in state space (cos/sin/w)
            F = figure(65);
            clf
            hold on
            jump_pt = [];
            for i = 1:length(obj.osm)
                osn = obj.osm{i};
                for j = 1:length(osn.sim)
                    traj_curr = osn.sim{j};
                    
                    
                    if (i==1) && (j == 1)
                        plot3(traj_curr.xl(:, 1), traj_curr.xl(:, 2), traj_curr.xl(:, 3), 'c', 'DisplayName', 'Trajectories')  
%                         scatter3(traj_curr.xl(end, 1), traj_curr.xl(end, 2), traj_curr.xl(end, 3), 30, 'ok', 'DisplayName', 'State Transition')  
                        jump_pt = [jump_pt; traj_curr.xl(end, :)];
                    else
                        plot3(traj_curr.xl(:, 1), traj_curr.xl(:, 2), traj_curr.xl(:, 3), 'c', 'HandleVisibility', 'Off')  
%                         scatter3(traj_curr.xl(end, 1), traj_curr.xl(end, 2), traj_curr.xl(end, 3), 30, 'ok','HandleVisibility', 'Off')  
                        jump_pt = [jump_pt; traj_curr.xl(end, :)];
%                         scatter3(traj_curr.xl(1, 1), traj_curr.xl(1, 2), traj_curr.xl(1, 3), 100, 'sk')
                    end
                    %symmetry of system dynamics (c,s,w) <-> (c,-s,-w)
                    plot3(traj_curr.xl(:, 1), -traj_curr.xl(:, 2), -traj_curr.xl(:, 3), 'c', 'HandleVisibility', 'Off')  
%                     scatter3(traj_curr.xl(end, 1), -traj_curr.xl(end, 2), -traj_curr.xl(end, 3), 30, 'ok','HandleVisibility', 'Off')  
                    jump_pt = [jump_pt; [traj_curr.xl(end, 1), -traj_curr.xl(end, 2:3)]];
                end
            end
            scatter3(jump_pt(:, 1),jump_pt(:, 2),jump_pt(:, 3), 30, 'ok', 'DisplayName', 'State Transition')  
            scatter3(-1, 0, 0, 400, 'rs', 'LineWidth', 2, 'DisplayName', 'Setpoint')
            
            %plot the circles
            if nargin > 1
                theta = linspace(0, 2*pi, 300);
                circ_x = cos(theta);
                circ_y = sin(theta);
                circ_z = zeros(size(theta));
                plot3(circ_x, circ_y, w_range_0(1) + circ_z, 'k', 'DisplayName', 'Initial Condition')
                plot3(circ_x, circ_y, w_range_0(2) + circ_z, 'k', 'HandleVisibility', 'Off')
                plot3(circ_x, circ_y, -w_max_est+circ_z, 'r--', 'DisplayName', 'Upper Bound')
                plot3(circ_x, circ_y, +w_max_est+circ_z, 'r--', 'HandleVisibility', 'Off')
            end
            
        %     xlabel('time', 'FontSize', FS_axis);
            xlabel('$cos(\theta)$', 'interpreter', 'latex', 'FontSize', obj.FS_axis);
            ylabel('$sin(\theta)$', 'interpreter', 'latex', 'FontSize', obj.FS_axis);
            zlabel('$\omega$', 'interpreter', 'latex', 'FontSize', obj.FS_axis);
            title('Pendulum Swing-up State Space', 'FontSize',  obj.FS_title);
            view(3)
            
            legend('location', 'northwest')
        end
        
        function F = cost_plot(obj, w_max_est)
            %COST_PLOT plot objective (w^2) of  trajectories in time
            F = figure(66);
            clf
            hold on
            for i = 1:length(obj.osm)
                osn = obj.osm{i};
                for j = 1:length(osn.sim)
                    traj_curr = osn.sim{j};
                    plot(traj_curr.t, traj_curr.x(:, 2).^2, 'c')
                    scatter(traj_curr.t(end), traj_curr.x(end, 2).^2, 30, 'ok')
                end
            end

            plot(xlim, w_max_est^2*[1;1], 'r--', 'LineWidth', 2)
            
        %     xlabel('time', 'FontSize', FS_axis);
            xlabel('$t$', 'interpreter', 'latex', 'FontSize', obj.FS_axis);
            ylabel('$\omega^2$', 'interpreter', 'latex', 'FontSize', obj.FS_axis);
%             zlabel('$\omega$', 'interpreter', 'latex', 'FontSize', obj.FS_axis);
            title('Angular Velocity Squared', 'FontSize',  obj.FS_title);
%             view(3)
        end
        
        function F = aux_plot(obj)
            %AUX_PLOT plot auxiliary function (v_i) of  trajectories in time
            F = figure(67);
            clf
            hold on
            jump_pt = [];
            for i = 1:length(obj.osm)
                osn = obj.osm{i};
                for j = 1:length(osn.sim)
                    traj_curr = osn.sim{j};
                    plot(traj_curr.t, traj_curr.v, 'c')
                    if j < length(osn.sim)
%                         scatter(traj_curr.t(end), traj_curr.v(end), 30, 'ok')
                        jump_pt = [jump_pt; [traj_curr.t(end), traj_curr.v(end)]];
                    end
                end
            end
            
            scatter(jump_pt(:, 1), jump_pt(:, 2), 30, 'ok')

%             plot(xlim, w_max_est^2*[1;1], 'r--', 'LineWidth', 2)
            
        %     xlabel('time', 'FontSize', FS_axis);
            xlabel('$t$', 'interpreter', 'latex', 'FontSize', obj.FS_axis);
            ylabel('$v(t)$', 'interpreter', 'latex', 'FontSize', obj.FS_axis);
%             zlabel('$\omega$', 'interpreter', 'latex', 'FontSize', obj.FS_axis);
            title('Auxiliary Function', 'FontSize',  obj.FS_title);
%             view(3)
        end
        
        
    end
end

