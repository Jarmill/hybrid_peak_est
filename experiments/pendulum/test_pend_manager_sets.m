
% PM.k = 2;
rng(22, 'twister')

% w_range = [-0.75; 0.75];
w_max_0 = 0.75;
w_range = [1; -1]*w_max_0;


SOLVE = 1;
SAMPLE = 1;
PLOT_NONNEG = 1;
PLOT = 0;
if SOLVE
%there was an optimal setting. What was it?
%maybe that was with starting at pi
%     Npartitions = 3;
%     order = 2;
%     %w_max_est = 2.050828859759454

    Npartitions = 0;
%     order = 1;
    %w_max_est =[2.11514442243362]
    
%     order = 2;
    %w_max_est = 2.05121731019132

    order = 3;
%     w_max_est = [2.04940053118856];
    
%     order = 4;
    %w_max_est = [2.04909496580124];
    
%     Npartitions = 4;
%     order = 3;
    % [2.04918210133319]

    
%     Npartitions = 6;
%     order = 1;
    
    PM = pend_manager();
    supp0 = [PM.vars.x(3)^2 <= w_max_0^2];
    PM = PM.make_manager(supp0, Npartitions);

    [PM, w_max_est, sol] = PM.solve(order);
    [op, mom_out, corner_out] = PM.PM.recover();
    op_all = all(op);
    

    %2 partitions at order 2
    %w_max_est = 1.901106516822999
end
if SAMPLE
    ps_func = @() PM.pend_sample_circ(w_range);
%     Nsample = 150;
    Nsample = 200;
%     Nsample = 20;
    ps = struct('N', Nsample, 'init', ps_func);
    
%     osn = PM.pend_sim_nonneg(x0, 10);
    
    osm = PM.pend_sample_multi(ps, 10);
end
if PLOT_NONNEG
    
    PP = pend_plotter(osm);
%     PP.state_plot_trig(w_range, w_max_est);
%     PP.nonneg_jump();
%     PP.nonneg_loc();
    PP.cost_plot(w_max_est);
    PP.track_plot_trig();


    
%     FS_title = 16;
%     FS_axis = 16;
%     
%     %jump nonnegativity
%     figure(60)
%     clf
%     hold on
%     title('Jump Nonnegativity', 'FontSize', FS_title)
%     ylabel("$v_i(x) - v(R_{i \rightarrow i'} x)$", 'interpreter', 'latex', 'FontSize', FS_axis);
%     xlabel('time', 'FontSize', FS_axis)
%     for i = 1:length(osm)
%         osn = osm{i};
%         for j = 1:length(osn.jump)
%             j_curr = osn.jump{j};            
%             stem(j_curr.t, j_curr.nonneg, 'c') 
%         end
%     end
%     
%     %location nonnegativity
%     figure(61)
%     clf
%     nonneg_title = {'Initial Value', 'Decrease in Value', 'Cost Proxy'};
%     ax_loc = {'$\gamma - v(x)$', '$-Lfv(x)$', '$v(x) - p(x)$'};
%     for k = 1:3
%         subplot(3, 1, k)
%         hold on
%         xlabel('time', 'FontSize', FS_axis)
%         ylabel(ax_loc{k}, 'interpreter', 'latex', 'FontSize', FS_axis);
%         title(nonneg_title{k}, 'FontSize', FS_title);   
%     end
%     for i = 1:length(osm)
%             osn = osm{i};
%         for j = 1:length(osn.sim)
%             traj_curr = osn.sim{j};
%             for k = 1:3            
%                 subplot(3, 1, k)
%                 plot(traj_curr.t, traj_curr.nonneg(:, k), 'c')
%                 scatter(traj_curr.t(1), traj_curr.nonneg(1, k), 100, 'ok')
%                 scatter(traj_curr.t(end), traj_curr.nonneg(end, k), 100, 'ok')
%             end
%         end  
%     end
%     
%     %Time Plot
%     figure(62)
%     clf
%     hold on
%     for i = 1:length(osm)
%         osn = osm{i};
%         for j = 1:length(osn.sim)
%             traj_curr = osn.sim{j};
%             plot3(traj_curr.t, traj_curr.xl(:, 2), -traj_curr.xl(:, 1), 'c')    
%             scatter3(traj_curr.t(end), traj_curr.xl(end, 2), -traj_curr.xl(end, 1), 100, 'ok')
%         end
%     end
%     xlabel('time', 'FontSize', FS_axis);
%     ylabel('$sin(\theta)$', 'interpreter', 'latex', 'FontSize', FS_axis);
%     zlabel('$-cos(\theta)$', 'interpreter', 'latex', 'FontSize', FS_axis);
%     title('Pendulum Swing-up Simulation', 'FontSize',  FS_title);
%     view(3)
%         
%     %State/Energy Plot
%     figure(63)
%     clf
%     for i = 1:length(osm)
%         osn = osm{i};
%         for j = 1:length(osn.sim)
%             traj_curr = osn.sim{j};
%             %theta
%             subplot(3, 1, 1)
%             hold on
%             theta_wrap = wrapTo2Pi(traj_curr.x(:, 1));
%             plot(traj_curr.t, theta_wrap, 'c')
%             scatter(traj_curr.t(1), theta_wrap(1), 40, 'ok')
%             scatter(traj_curr.t(end), theta_wrap(end), 40, 'ok')
%             
% 
%             %w
%             subplot(3, 1, 2)
%             hold on
%             plot(traj_curr.t, traj_curr.x(:, 2), 'c')
%             scatter(traj_curr.t(1), traj_curr.x(1, 2), 40, 'ok')
%             scatter(traj_curr.t(end), traj_curr.x(end, 2), 40, 'ok')
%             
% 
%             %energy
%             subplot(3, 1, 3)
%             hold on
%             plot(traj_curr.t, traj_curr.E, 'c')
%             scatter(traj_curr.t(1), traj_curr.E(1), 40, 'ok')
%             scatter(traj_curr.t(end), traj_curr.E(end), 40, 'ok')
% 
%         end
%     end
%     %cleanup
%     subplot(3, 1, 1)
%     ylim([0, 2*pi])
%     plot(xlim, [pi, pi], ':k')
%     title('Angle', 'FontSize', FS_title)
%     
%     
%     subplot(3, 1, 2)
%     plot(xlim, [0, 0], ':k')
%     title('Angular Velocity', 'FontSize', FS_title)                    
%     xlabel('time')
%     
%     subplot(3, 1, 3)
%     plot(xlim, [1, 1], ':k')
%     xlabel('time')
%     title('Energy', 'FontSize', FS_title)
%     
%     
%     %Trig State Space Plot        
%     figure(64)
%     clf
%     hold on
%     for i = 1:length(osm)
%         osn = osm{i};
%         for j = 1:length(osn.sim)
%             traj_curr = osn.sim{j};
%             plot3(traj_curr.xl(:, 1), traj_curr.xl(:, 2), traj_curr.xl(:, 3), 'c')  
%             if j == 1
%                 scatter3(traj_curr.xl(1, 1), traj_curr.xl(1, 2), traj_curr.xl(1, 3), 100, 'ok')
%             end
%             scatter3(traj_curr.xl(end, 1), traj_curr.xl(end, 2), traj_curr.xl(end, 3), 50, 'ok')
%         end
%     end
%     scatter3(-1, 0, 0, 400, 'rs', 'LineWidth', 2)
% %     xlabel('time', 'FontSize', FS_axis);
%     xlabel('$cos(\theta)$', 'interpreter', 'latex', 'FontSize', FS_axis);
%     ylabel('$sin(\theta)$', 'interpreter', 'latex', 'FontSize', FS_axis);
%     zlabel('$\omega$', 'interpreter', 'latex', 'FontSize', FS_axis);
%     title('Pendulum Swing-up State Space', 'FontSize',  FS_title);
%     view(3)

end 
if PLOT
    os = PM.pend_sim(x0, 20);
    w_max_true = max(abs(os.theta(:, 2)));
    figure(2)
    clf
    subplot(3, 1, 1)
    plot(os.t, wrapTo2Pi(os.theta(:, 1)))
    ylim([0, 2*pi])
    hold on
    plot(xlim, [pi, pi], ':k')
    title('Angle')
    subplot(3, 1, 2)
    plot(os.t, os.theta(:, 2))
    title('Angular Velocity')
    hold on
    plot(xlim, [0, 0], ':k')
    subplot(3, 1, 3)
    E = PM.energy(os.x')';
    plot(os.t, E)
    title('Energy')
    hold on
    plot(xlim, [1, 1], ':k')
    xlabel('time')
end