% the prajna disturb_cube example (Ex. 2 from
% https://viterbi-web.usc.edu/~jdeshmuk/teaching/cs699-fm-for-cps/Papers/A5.pdf)
% is extremely ill-conditioned when using the monomial basis.
%
% This script aims to find an alternative set of parameters that are better
% conditioned to yield demonstrative and feasible solutions for hybrid peak
% estimation
%
% There is stochastic time-dependent noise in this example. This script investigates whether the
% auxiliary function v(t, x) decreases along noise-free trajectories



SETUP = 1;
SOLVE = 0;
SAMPLE = 1;
PLOT = 0;


%% setup the locations and guards
if SETUP
    
    mset clear
    mpol('t');
    mpol('x', 3, 1);

    vars = struct('t', t, 'x', x);

%     Tmax = 5;
%     Tmax = 10;
%     Tmax = 15;
    Tmax = 20;
%     Tmax = 0;

    %box 1 (no control)
    r2 = sum(x.^2);
    
    K = 0.25;
    critfw = K*x(1)^2 + x(2)^2 + x(3)^2;
    
    R0 = 0.2;
%     R1 = 1;
    R1 = 1.5;
    
    sigma = 0.1;
    
    X01 = [r2 == R0^2];
%     X01 = [R0;0;0];
    
    L = R1*[1;1;1];
        
    lsupp1 = chance_support(vars);
    lsupp1.bound_type = 'mean';
    lsupp1 = lsupp1.set_box(L);
    lsupp1.X = [lsupp1.X; critfw <= R1^2];
    lsupp1.X_init = X01;
    lsupp1.Tmax = Tmax;
%     lsupp1.disturb = w^2<=1;
    
%     p1 = x(1)^2;
    p1 = [];
%     f1 = [x(2); -x(1) + x(3); x(1) + (2*x(2) + 3*x(3))*(1+x(3)^2) + w];
f1 = [x(2); (-x(1) + x(3)); x(1) + (2*x(2) + 3*x(3))*(1+x(3)^2)];
g1 = [0; 0; sigma];
%     f1 = [x(2); -x(1) + x(3); x(1) + (2*x(2) + 3*x(3)) + w];

    dyn1 = struct('f', f1, 'g', g1);

    loc1 = location_sde(lsupp1, dyn1, p1, 1);
    
        
    %box 2 (control)
    X02 = [];
    
    lsupp2 = chance_support(vars);
    lsupp2.bound_type = 'mean';
    lsupp2 = lsupp2.set_box(L);
    lsupp2.X = [lsupp2.X; r2 >= R0^2];
    lsupp2.X_init = X02;
    lsupp2.Tmax = Tmax;
%     lsupp2.disturb = w^2<=1;
    
    
    % X02 = [t == 0;  (x - [-0.4; 0]).^2 <= 0.01];
    
%     slow2 = 0.75;
    slow2 = 1;
    f2 = [x(2); (-x(1) + x(3)); -x(1) - (2*x(2) + 3*x(3))];
    g2 = [0; 0; sigma];
    p2 = x(1)^2;
%     p2 = [x(1)^2; x(2)^2];

    dyn2 = struct('f', f2, 'g', g2);

    loc2 = location_sde(lsupp2, dyn2, p2, 2);
    
    %guards
    R  = x; %reset map
    
    
    %possible bug: time not defined in support in the reset maps
    Xfw = [critfw == R1^2; t*(1-t)>=0];
    
    critbk = r2;
    Xbk = [critbk == R0^2; t*(1-t)>=0];   
    
    gfw = guard(1, vars, loc1, loc2, Xfw, R);
    gbk = guard(2, vars, loc2, loc1, Xbk, R);
    
    Zeno_cap = Inf;
    gfw.zeno_cap = Zeno_cap;
    gbk.zeno_cap = Zeno_cap;
end

%% solve the system and get peak estimates
if SOLVE
    
    PM =  peak_manager_hy({loc1, loc2}, {gfw, gbk});

%SDE with sigma=0.1
% order = 1; % 2.2500
% order = 2; % 0.6904 
% order = 3; %0.5174 
% order = 4; %0.4691
order = 5; %0.4525 
% order = 4;

 p_order = zeros(length(orderlist), 1);
    time_order = zeros(length(orderlist), 1);
     for i = 1:length(orderlist)
    
    %     [objective, mom_con, supp_con, len_dual] =  PM.cons(order);
      % tic
            [sol, PM] = PM.run(orderlist(i), Tmax);
            time_order(i) = sol.solver_time;
            p_order(i) = sol.obj_rec;
    %     [sol, = PM.run(order) ;
    %     fprintf('abs(x1) bound: %0.4f \n', sqrt(sol.obj_rec))
    fprintf('bound: %0.4f \n', (sol.obj_rec))
        p_est = sqrt(sol.obj_rec);
        % [rr, mm, cc] = PM.recover();
    
        obj_rec = sol.obj_rec;
        save('stoch_cube.mat', 'time_order', 'p_order', 'orderlist');
     end
    % [rr, mm, cc] = PM.recover();
end

%% sample trajectories
if SAMPLE
    rng(24, 'twister')
%     smp1 = struct('x', @() sphere_sample(1, 3)'*R0);
% smp1 = struct('x', @() sphere_sample(1, 3)'*(1e-3+R0));

Nbatch = 1;
Nsample = 1000;
x0_list = zeros(3, Nbatch);
    osm_list = cell(Nbatch, 1);
    osd_list = cell(Nbatch, 1);

    for i = 1:Nbatch
        % x0 = sphere_sample(1, 3)'*R0;
        x0 = [0.5; 0.5; 0.5];
        x0_list(:, i) = x0;
        smp1 = struct('x', @() x0);
    
        smp2 = struct('x', []);
        
        
        LS1 = sampler_sde_uncertain(loc1, smp1);
        
        LS2 = sampler_sde_uncertain(loc2, smp2);
        
        
        LS1.mu = 0.01;
        LS2.mu = 0.01;
    
    %     LS1.mu = 0.03;
    %     LS2.mu = 0.03;
        
        HS = sampler_hy({LS1, LS2}, {gfw, gbk});
        
        %     osh = HS.sample_traj(0, [0;0;0.03], 1, 5);
    %     Nsample = ;
    %     Nsample = 20;
        
        % Nsample = 5000;
        % Nsample = 5;
    
        
    
    
        [osm_list{i}, osd_list{i}] = HS.sample_traj_multi(Nsample, 1);
        
        

        % t_end = cellfun(@(o) o.t_end, osm);
    end

    % osd = struct('locations', cell(2, 1));
    osd = struct;
    osd.locations = cell(2, 1);    
    osd.locations{1} = osd_list{1}.locations{1};
    osd.locations{2} = osd_list{1}.locations{2};
    for i = 2:Nbatch
        osd.locations{1} = [osd.locations{1}; osd_list{i}.locations{1}];
        osd.locations{2} = [osd.locations{2}; osd_list{i}.locations{2}];
    end
    
%     osm = SMP.sample_traj_multi(Ntraj, lsupp.Tmax);
    % save('traj/stoch_cube_traj_big.mat', 'osm', 'osd')
    save('traj/stoch_cube_traj_big.mat', 'osm_list', 'osd_list')
else
    load('traj/stoch_cube_traj_big.mat');
    Ntraj = length(osm);
end
    

%% plot trajectories
if PLOT
    CPlot = cube_plotter(osm, osd, R0, R1);
    CPlot.cube_plot();
%     CPlot.cube_plot(sqrt(sol.obj_rec));
%     CPlot.nonneg_loc();
%     CPlot.aux_plot();
%     CPlot.nonneg_jump();

%     CPlot.nonneg_jump();
%     CPlot.objective_plot();
end