rng(3, 'twister');
param = struct('R0', 0.2, 'R1', 1.5, 'K', 0.25, 'sigma', 0.1, 'dt', 1e-3, 'L', 1.5);

% x0 = [0.5; 0.5; 0.5];

Tmax = 5;

%% sample
% osd = stoch_cube_sampler_single(x0, Tmax, param);



Nbatch = 300;
osd_list = cell(Nbatch, 1);
max_mean = zeros(Nbatch, 1);
t_raw = 0:1e-3:Tmax;
for j = 1:Nbatch
    x_2_trace = zeros(k,length(t_raw));
    Ntraj = 100;
    osd.locations = {{}, {}};
    osd.guards = {{}, {}};
    % x0 = [0; 0; 0];
    x0 = sphere_sample(1, 3)'*param.R0;
    for k = 1:Ntraj
        osd_curr = stoch_cube_sampler_single(x0, Tmax, param);
        
        for i = 1:length(osd.locations)
            osd.locations{i} = vertcat(osd.locations{i}, osd_curr.locations{i});
        end
        
        for i = 1:length(osd.guards)
            osd.guards{i} = vertcat(osd.guards{i}, osd_curr.guards{i});
        end
    
        %analyze trajectory
        for i = 1:length(osd_curr.locations{2})
            curr = osd_curr.locations{2}{i};
            x_2_trace(k, 1 + uint64(1e3*(curr.t))) = curr.x(2:end, 1)';
        end        
    end
    p_trace = x_2_trace.^2;
    mean_trace = sum(p_trace, 1)/Ntraj;
    max_mean(j) = max(mean_trace);
    osd_list{j} = osd;
end


save('traj_stoch_cube_big.mat', 'osd_list', 'param', 'max_mean');

%% plot
    CPlot = cube_plotter([], osd, param.R0, param.R1);
    CPlot.cube_plot(sqrt(0.4525));

