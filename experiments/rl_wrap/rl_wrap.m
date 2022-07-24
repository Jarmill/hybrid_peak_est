% A single-location case where the left side jumps to the bottom, and the
% right side jumps to the top. Both the top and bottom are purely
% inward-facing.
%
% There is no noise in this example. This script investigates whether the
% auxiliary function v(t, x) decreases along noise-free trajectories



SETUP = 1;
SOLVE = 0;
SAMPLE = 1;
PLOT = 1;


%% setup the locations and guards
if SETUP
    
    mset clear
    mpol('t');
    mpol('x', 2, 1);
%     mpol('w', 1, 1);
    w=0;
    vars = struct('t', t, 'x', x);

%     Tmax = 5;
%     Tmax = 10;
    Tmax = 15;
%     Tmax = 20;
%     Tmax = 0;
% Tmax = 3;

    %box 1 (no control)
    
    R0 = 0.2;
%     C0 = [-0.75; 0.75];
    C0 = [0.75; 0.75];
    R1 = 1;
%     R1 = 1.5;
%     R1=2;
    
    X01 = [sum((x-C0).^2) <= R0^2];
%     X01 = [R0;0;0];
    
    L = R1*[1;1];
        
    lsupp1 = loc_support(vars);
    lsupp1 = lsupp1.set_box(L);
%     lsupp1.X = [lsupp1.X];
    lsupp1.X_init = X01;
    lsupp1.Tmax = Tmax;
%     lsupp1.disturb = w^2<=1;
    
    p1 = -x'*x;
%     p1 = [];

%     f1 = [-x(2) + x(1)*x(2)/2;
%          -x(2) - x(1) + x(1)^3/4];

    f1 = [-x(2) + x(1)*x(2) + 0.5;
         -x(2) - x(1) + x(1)^3];

    loc1 = rl_location(lsupp1, {f1}, p1, 1);

    
        
  %guards

    %reset maps
    Rleft = [-x(2); x(1)];     %left->bottom
    Rright = [x(2); x(1)]; %right->top
    
    
    %possible bug: time not defined in support in the reset maps
    Xleft  = [x(2)^2 <= R1^2; x(1)==-R1; t*(1-t)>=0];   
    Xright = [x(2)^2 <= R1^2; x(1)==R1; t*(1-t)>=0];   
    
    gleft  = guard(1, vars, loc1, loc1, Xleft, Rleft);
    gright = guard(2, vars, loc1, loc1, Xright, Rright);
    
    
    Zeno_cap = 5;
    gleft.zeno_cap = Zeno_cap;
    gright.zeno_cap = Zeno_cap;
end

%% solve the system and get peak estimates
if SOLVE
    
    PM =  peak_manager_hy({loc1}, {gleft, gright});

    %T = 5
    order = 1; %0
    order = 2; %0
%     order = 3; % 0
    order = 4; % 0
    order=5; % -0.010209422355083
    order = 6; % -0.019157038491643

%     [objective, mom_con, supp_con] =  PM.cons(order);
    [sol, PM] = PM.run(order, Tmax);    
%     sol = PM.run(order) ;
    obj_rec = sol.obj_rec;
%     fprintf('abs(x1) bound: %0.4f \n', sqrt(sol.obj_rec))
fprintf('-norm(x)^2 bound: %0.4f \n', (sol.obj_rec))
    p_est = sqrt(sol.obj_rec);
    [rr, mm, cc] = PM.recover();
end

%% sample trajectories
if SAMPLE
    rng(26, 'twister')
        smp1 = struct('x', @() sphere_sample(1, 2)'*R0 + C0);
%     smp1 = struct('x', @() C0);

%     smp2 = struct('x', []);
    
    
    LS1 = sampler_uncertain(loc1, smp1);
    
%     LS2 = sampler_uncertain(loc2, smp2);
    
    
    LS1.mu = 0.03;
%     LS2.mu = 0.03;
    
    HS = sampler_hy({LS1}, {gleft, gright});
    
    %     osh = HS.sample_traj(0, [0;0;0.03], 1, 5);
%     Nsample = 1;
% Nsample = 10;
%     Nsample = 20;
%     Nsample = 50;
    Nsample = 100;
%     Nsample = 5;


    [osm, osd] = HS.sample_traj_multi(Nsample, Tmax);
    
    t_end = cellfun(@(o) o.t_end, osm);
    
end

%% plot trajectories
if PLOT

    RPlot = rl_plotter(osm, osd, R0, C0);
    RPlot.rl_plot();
    RPlot.rl_plot(-0.019157038491643);
%     figure(3);
%     osc = osd.locations{1};
%     for j = 1:length(osc)
%         hold on
%         plot(osc{j}.x(:, 1), osc{j}.x(:, 2), 'c')
%     end
%     CPlot = cube_plotter(osm, osd, R0, R1);
%     CPlot.cube_plot(sqrt(sol.obj_rec));
%     CPlot.nonneg_loc();
%     CPlot.aux_plot();
%     CPlot.nonneg_jump();
% 
% %     CPlot.nonneg_jump();
% %     CPlot.objective_plot();
end