
% PM.k = 2;

% x0 = [pi; -1.5];
% x0 = [pi; 0.5];
% x0 = [pi/2; 0];
% x0 = [1e-2; 0];
% x0 = [pi/2; 0]; %optimal
% x0 = [pi/4; 0]; %maybe not optimal
% x0 = [pi/8; 0]; %maybe not optimal
% x0 = [pi/16; 0]; %maybe not optimal
% x0 = [pi/32; 0]; %maybe not optimal
% x0 = [1e-2; 0]; %maybe not optimal
SOLVE = 0;
SAMPLE = 1;
PLOT_NONNEG = 1;
PLOT = 1;
if SOLVE
%there was an optimal setting. What was it?
%maybe that was with starting at pi
    
%%Npartitions = 0

%     Npartitions = 0;
%     order = 4;

%     %w_max_est = 2.046643896001510    
%     Npartitions = 0;
%     order = 4;
    
%     %w_max_est = 1.960179247920592  
%     Npartitions = 0;
%     order = 5;

%%Npartitions = 1    
%     %wmax_est = 1.990651517936928
%     Npartitions = 1;
%     order = 3;
    
%     %wmax_est = 1.956749536836139
%     %with the bugged guard layout this was optimal. what happened?
%     %approximately optimal (compare to [1.95125300492727])
%     %measures are approx rank-1 (op = [1;1;1;1;1;1;1;1])
%     Npartitions = 1;
%     order = 4;

%Npartitions = 2
    %wmax_est = [1.949358851405596]
%     Npartitions = 2;
%     order = 1;


    Npartitions = 3;
    order = 2;




%     %wmax_est = [1.787994405288204 invalid]
%     Npartitions = 2;
%     order = 3;

% % Npartitions = 3
% %     wmax_est = [2.04939008902368]
%     Npartitions = 3;
%     order = 1;

%     Npartitions = 3;
%     order = 2;

% %     wmax_est = [1.79728275516481]
%     Npartitions = 3;
%     order = 2;

% %     wmax_est = [1.63333845511395]
%     Npartitions = 3;
%     order = 3;


% % Npartitions = 4
% %     wmax_est = [2.04939013063602] 
%     Npartitions = 4;
%     order = 1;

% % %     wmax_est = [1.54782484972908] 
%     Npartitions = 4;
%     order = 2;

% % Npartitions = 5
% %     wmax_est = [2.04939013637590] 
%     Npartitions = 5;
%     order = 1;


Tmax = 10;
    


w_range = 0.5*[-1; 1];

    PM = pend_manager();
    supp0 = x0;
%     supp0 = [PM.vars.t == 0; PM.vars.x == PM.trig_lift(x0)];
    supp0 = [PM.vars.x(3)^2 <= w_range(2)^2];
    PM = PM.make_manager(supp0, Npartitions, Tmax);

    [PM, w_max_est, sol] = PM.run(order);
%     [sol, PM] = PM.run(order);
    [op, mom_out, corner_out] = PM.PM.recover();
    op_all = all(op);
    %without partitions at order 2:
    %w_max_est = 2.049354326410954

    %2 partitions at order 2
    %w_max_est = 1.901106516822999
end
if SAMPLE
        
%     osn = PM.pend_sim_nonneg(x0, 10);



    Nsample = 10;
    ps_func = @() PM.pend_sample_point_box(c_range, w_range);
    ps = struct('N', 10, 'init', ps_func);
    osn = PM.pend_sample_multi(ps, Tmax);

end

if PLOT_NONNEG
    FS_title = 16;
    FS_axis = 16;
    figure(60)
    clf
    hold on
    title('Jump Nonnegativity', 'FontSize', FS_title)
    ylabel('$v(x) - v(Rx)$', 'interpreter', 'latex', 'FontSize', FS_axis);
    xlabel('time', 'FontSize', FS_axis)

    for i = 1:length(osn)
    for j = 1:length(osn{i}.jump)
        j_curr = osn{i}.jump{j};            
        stem(j_curr.t, j_curr.nonneg, 'c') 
    end
    end
    hold off

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
%     %     title(['Loc ', num2str(i), ' ', nonneg_title{k}], 'FontSize', FS_title)
% 
%     end
% 
%     for j = 1:length(osn.sim)
%         traj_curr = osn.sim{j};
%         for k = 1:3            
%             subplot(3, 1, k)
%             plot(traj_curr.t, traj_curr.nonneg(:, k), 'c')
%             scatter(traj_curr.t(1), traj_curr.nonneg(1, k), 100, 'ok')
%             scatter(traj_curr.t(end), traj_curr.nonneg(end, k), 100, 'ok')
%         end
%     end  

    figure(62)
    clf
    hold on
    for i = 1:length(osn)
    for j = 1:length(osn{i}.sim)
        traj_curr = osn{i}.sim{j};
        plot3(traj_curr.t, traj_curr.xl(:, 2), -traj_curr.xl(:, 1), 'c')    
        scatter3(traj_curr.t(end), traj_curr.xl(end, 2), -traj_curr.xl(end, 1), 100, 'ok')
    end
    end
    xlabel('time', 'FontSize', FS_axis);
    ylabel('$cos(\theta)$', 'interpreter', 'latex', 'FontSize', FS_axis);
    zlabel('$cos(\theta)$', 'interpreter', 'latex', 'FontSize', FS_axis);
    title('Pendulum Swing-up Simulation', 'FontSize',  FS_title);
    view(3)

end 
if PLOT
%     os = PM.pend_sim(x0, 20);
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