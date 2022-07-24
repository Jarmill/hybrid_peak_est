%after running deterministic_cube.m

load("traj_deterministic_cube.mat")

%locations
for i = 1:length(PM.loc)
    loc_curr = PM.loc{i};
    for j = 1:length(osd.locations{i})
        osc = osd.locations{i}{j};
        osd.locations{i}{j}.nonneg = PM.loc{i}.nonneg(osc.t, osc.x, [], []);
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


%% 
    CPlot = cube_plotter(osm, osd, R0, R1);
%     CPlot.cube_plot(sqrt(sol.obj_rec));
    CPlot.nonneg_loc();
    CPlot.aux_plot();
    CPlot.nonneg_jump();
