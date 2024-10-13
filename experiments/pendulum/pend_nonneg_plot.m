FS_axis = 14;
FS_title = 14;

osn = PM.pend_sim_nonneg(x0, 10);

figure(60)
clf
hold on
title('Jump Nonnegativity', 'FontSize', FS_title)
ylabel('$v(x) - v(Rx)$', 'interpreter', 'latex', 'FontSize', FS_axis);
xlabel('time', 'FontSize', FS_axis)

for j = 1:length(osn.jump)
    j_curr = osn.jump{j};            
    stem(j_curr.t, j_curr.nonneg, 'c') 
end
hold off

figure(61)
clf
nonneg_title = {'Initial Value', 'Decrease in Value', 'Cost Proxy'};
ax_loc = {'$\gamma - v(x)$', '$-Lfv(x)$', '$v(x) - p(x)$'};
for k = 1:3
    subplot(3, 1, k)
    hold on
    xlabel('time', 'FontSize', FS_axis)
    ylabel(ax_loc{k}, 'interpreter', 'latex', 'FontSize', FS_axis);
    title(nonneg_title{k}, 'FontSize', FS_title);
%     title(['Loc ', num2str(i), ' ', nonneg_title{k}], 'FontSize', FS_title)
    
end

for j = 1:length(osn.sim)
    traj_curr = osn.sim{j};
    for k = 1:3            
        subplot(3, 1, k)
        plot(traj_curr.t, traj_curr.nonneg(:, k), 'c')
        scatter(traj_curr.t(1), traj_curr.nonneg(1, k), 100, 'ok')
        scatter(traj_curr.t(end), traj_curr.nonneg(end, k), 100, 'ok')
    end
end  

figure(62)
clf
hold on
for j = 1:length(osn.sim)
    traj_curr = osn.sim{j};
    plot3(traj_curr.t, traj_curr.xl(:, 2), -traj_curr.xl(:, 1), 'c')    
    scatter3(traj_curr.t(end), traj_curr.xl(end, 2), -traj_curr.xl(end, 1), 100, 'ok')
end
xlabel('time', 'FontSize', FS_axis);
ylabel('$cos(\theta)$', 'interpreter', 'latex', 'FontSize', FS_axis);
zlabel('$cos(\theta)$', 'interpreter', 'latex', 'FontSize', FS_axis);
title('Pendulum Swing-up Simulation', 'FontSize',  FS_title);
view(3)

