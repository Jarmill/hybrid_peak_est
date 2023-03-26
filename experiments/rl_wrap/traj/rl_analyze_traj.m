%after running deterministic_cube.m

% load("rl_traj_15.mat")
% load("rl_traj_3.mat")

%locations
for i = 1:length(PM.loc)
    loc_curr = PM.loc{i};
    for j = 1:length(osd.locations{i})
        osc = osd.locations{i}{j};
%         osd.locations{i}{j}.nonneg = PM.loc{i}.nonneg(osc.t, osc.x, [], []);
        osd.locations{i}{j}.v = PM.loc{i}.v_eval(osc.t, osc.x, []);
        osd.locations{i}{j}.obj = PM.loc{i}.obj_eval(osc.t, osc.x);
    end
end

for i = 1:length(PM.guards)
    g_curr = PM.guards{i};
    for j = 1:length(osd.guards{i})
        osg = osd.guards{i}{j};
        osd.guards{i}{j}.nonneg = g_curr.nonneg(osg.t, osg.x);
    end
end

 
%% output plots
    RPlot = rl_plotter(osm, osd, R0, C0, Cu);
    RPlot.rl_plot(obj_rec);
%     CPlot.cube_plot(sqrt(sol.obj_rec));
%     RPlot.nonneg_loc();
%     RPlot.aux_plot();
%     RPlot.nonneg_jump();
