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
PLOT = 1;


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
    
    Zeno_cap = 5;
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

    [objective, mom_con, supp_con] =  PM.cons(order);
%     [sol, PM] = PM.run(order, Tmax);    
    sol = PM.run(order) ;
%     fprintf('abs(x1) bound: %0.4f \n', sqrt(sol.obj_rec))
fprintf('x1^2 bound: %0.4f \n', (sol.obj_rec))
    p_est = sqrt(sol.obj_rec);
    [rr, mm, cc] = PM.recover();
end

%% sample trajectories
if SAMPLE
    rng(25, 'twister')
%     smp1 = struct('x', @() sphere_sample(1, 3)'*R0);
smp1 = struct('x', @() sphere_sample(1, 3)'*(1e-3+R0));
%     smp1 = struct('x', @() [R0; 0; 0]);

    smp2 = struct('x', []);
    
    
    LS1 = sampler_sde_uncertain(loc1, smp1);
    
    LS2 = sampler_sde_uncertain(loc2, smp2);
    
    
    LS1.mu = 0.01;
    LS2.mu = 0.01;

%     LS1.mu = 0.03;
%     LS2.mu = 0.03;
    
    HS = sampler_hy({LS1, LS2}, {gfw, gbk});
    
    %     osh = HS.sample_traj(0, [0;0;0.03], 1, 5);
    Nsample = 1;
%     Nsample = 20;
%     Nsample = 50;
%     Nsample = 5;


    [osm, osd] = HS.sample_traj_multi(Nsample, 1);
    
    t_end = cellfun(@(o) o.t_end, osm);
    
%     osm = SMP.sample_traj_multi(Ntraj, lsupp.Tmax);
    save('stoch_cube_traj.mat', 'osm', 'osd')
else
    load('stoch_cube_traj.mat');
    Ntraj = length(osm);
end
    

%% plot trajectories
if PLOT
    CPlot = cube_plotter(osm, osd, R0, R1);
%     CPlot.cube_plot();
    CPlot.cube_plot(sqrt(sol.obj_rec));
%     CPlot.nonneg_loc();
%     CPlot.aux_plot();
%     CPlot.nonneg_jump();

%     CPlot.nonneg_jump();
%     CPlot.objective_plot();
end