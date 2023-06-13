rng(3, 'twister');
param = struct('R0', 0.2, 'R1', 1.5, 'K', 0.25, 'sigma', 0.1, 'dt', 1e-3, 'L', 1.5);

% x0 = [0.5; 0.5; 0.5];
smp1 = struct('x', @() sphere_sample(1, 3)'*param.R0);
Tmax = 20;

%% sample
% osd = stoch_cube_sampler_single(x0, Tmax, param);

osd.locations = {{}, {}};
osd.guards = {{}, {}};
Ntraj = 200;
for k = 1:Ntraj
    x0 = smp1.x();
    osd_curr = stoch_cube_sampler_single(x0, Tmax, param);
    
    for i = 1:length(osd.locations)
        osd.locations{i} = vertcat(osd.locations{i}, osd_curr.locations{i});
    end
    
    for i = 1:length(osd.guards)
        osd.guards{i} = vertcat(osd.guards{i}, osd_curr.guards{i});
    end
    
    if mod(k, 10)==0
        disp(k)
    end
    
end

save('traj_stoch_cube_sphere.mat', 'osd', 'param');


%% plot
    CPlot = cube_plotter([], osd, param.R0, param.R1);
    CPlot.cube_plot();

