
% A single-location case where the left side jumps to the bottom, and the
% right side jumps to the top. Both the top and bottom are purely
% inward-facing.
%
% There is no noise in this example. This script investigates whether the
% auxiliary function v(t, x) decreases along noise-free trajectories



SETUP = 1;
SOLVE = 1;
SAMPLE = 0;
PLOT = 0;
% rl_analyze_traj;


%% setup the locations and guards
if SETUP
    
    mset clear
    mpol('t');
    mpol('x', 2, 1);
%     mpol('w', 1, 1);
    w=0;
    vars = struct('t', t, 'x', x);

%     Tmax = 4;
    Tmax = 5;
%     Tmax = 10;
%     Tmax = 15;
%     Tmax = 20;
%     Tmax = 0;
% Tmax = 3;

    %box 1 (no control)
    
    R0 = 0.2;
%     C0 = [-0.75; 0.75];
%     C0 = [0.75; -0.75];
    C0 = [0.5; -0.3];
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
    
%     p1 = -x'*x;

    Cu = [-0.5; -0.5];
    p1 = -sum((x-Cu).^2);

%     p1 = -sum(x);
%     p1 = [];

%     f1 = [-x(2) + x(1)*x(2)/2;
%          -x(2) - x(1) + x(1)^3/4];

    %original
    %k = 0.5;
    k = 0.6;
    f1 = [-x(2) + x(1)*x(2) + k;
         -x(2) - x(1) + x(1)^3];

    loc1 = rl_location(lsupp1, {f1}, p1, 1);

    
        
  %guards

    %reset maps
    Rleft = [1-x(2); x(1)];     %left->bottom
    Rright = [x(2); x(1)]; %right->top
    
    
    %possible bug: time not defined in support in the reset maps
    Xleft  = [x(2)^2 <= R1^2; x(1)==-R1; t*(1-t)>=0];   
    Xright = [x(2)^2 <= R1^2; x(1)==R1; t*(1-t)>=0];   
    
    gleft  = guard(1, vars, loc1, loc1, Xleft, Rleft);
    gright = guard(2, vars, loc1, loc1, Xright, Rright);
    
    
    Zeno_cap = 0;
%     Zeno_cap = 10;
%     Zeno_cap = 15;
    gleft.zeno_cap = Zeno_cap;
    gright.zeno_cap = Zeno_cap;
end

%% solve the system and get peak estimates
if SOLVE
    
    PM =  peak_manager_hy({loc1}, {gleft, gright});

%     %T = 5, N =5
%     order = 1; %0
%     order = 2; %0
% %     order = 3; % 0
%     order = 4; % 0
%     order=5; % -0.010209422355083
%     order = 6; % -0.019157038491643

%T = 5, N = 5
% order = 1; %0
%     order = 2; %0
%     order = 3; % -0.364417031631220
%     order = 4; % -0.525881280208672
%    order = 5; %-0.565907365740604
   % order = 6; %-0.572099964307251
orderlist = 1:6;
 % orderlist = 1;
 
    % orderlist = [1; 2; 1];
    % orderlist = 6;
    p_order = zeros(length(orderlist), 1);
    time_order = zeros(length(orderlist), 1);
     for i = 1:length(orderlist)
    
    %     [objective, mom_con, supp_con, len_dual] =  PM.cons(order);
      % tic
            [PM, sol] = PM.run(orderlist(i), Tmax);
            time_order(i) = sol.solver_time;
            p_order(i) = sol.obj_rec;
    %     [sol, = PM.run(order) ;
    %     fprintf('abs(x1) bound: %0.4f \n', sqrt(sol.obj_rec))
    fprintf('bound: %0.4f \n', (sol.obj_rec))
        p_est = sqrt(sol.obj_rec);
        % [rr, mm, cc] = PM.recover();
    
        obj_rec = sol.obj_rec;
        % save('rl_wrap.mat', 'time_order', 'p_order', 'orderlist');
     end
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
% Nsample = 30;
%     Nsample = 50;
    Nsample = 120;
%     Nsample = 5;


    [osm, osd] = HS.sample_traj_multi(Nsample, Tmax);
    
    t_end = cellfun(@(o) o.t_end, osm);
    
end

%% plot trajectories
if PLOT

    RPlot = rl_plotter(osm, osd, R0, C0, Cu);
    RPlot.rl_plot();
    
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