% the prajna disturb_cube example (Ex. 2 from
% https://viterbi-web.usc.edu/~jdeshmuk/teaching/cs699-fm-for-cps/Papers/A5.pdf)
% is extremely ill-conditioned when using the monomial basis.
%
% This script aims to find an alternative set of parameters that are better
% conditioned to yield demonstrative and feasible solutions for hybrid peak
% estimation
%
% There is no noise in this example. This script investigates whether the
% auxiliary function v(t, x) decreases along noise-free trajectories


%this performs a safety-margin analysis on the unsafe set:
% Ru = 0.25;
% Cu = [0.7; 0; 0.5];
% Xu = [Ru^2 - (x-Cu).^2; x(3) - Cu(3)]>=0;



SETUP = 1;
SOLVE = 1;
SAMPLE = 0;
PLOT = 0;


%% setup the locations and guards
if SETUP
    
    mset clear
    mpol('t');
    mpol('x', 3, 1);
%     mpol('w', 1, 1);
    w=0;
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
    
    X01 = [r2 == R0^2];
%     X01 = [R0;0;0];
    
    L = R1*[1;1;1];
        
    lsupp1 = loc_support(vars);
    lsupp1 = lsupp1.set_box(L);
    lsupp1.X = [lsupp1.X; critfw <= R1^2];
    lsupp1.X_init = X01;
    lsupp1.Tmax = Tmax;
%     lsupp1.disturb = w^2<=1;
    
%     p1 = x(1)^2;
%     p1 = [];

%     Ru = 0.25;
% Ru = 0.4;
% Ru = 0.3;
Ru = 0.25;
%     Cu = [0.7; 0; 0.5];
%     Cu = [0.55; 0; 0.5];
%     Cu = [-0.5; -0.3; 0.5];
% Cu = [-0.5; -0.5; 0.5];
% Cu = [-0.6; -0.6; 0.5]; %is safe at order 4
Cu = [-0.6; -0.5; 0.5]; %
p1 = [Ru^2 - sum((x-Cu).^2); x(3) - Cu(3)];

%     p1 = Ru^2 - (x-Cu).^2;
%     f1 = [x(2); -x(1) + x(3); x(1) + (2*x(2) + 3*x(3))*(1+x(3)^2) + w];
f1 = [x(2); (-x(1) + x(3)); x(1) + (2*x(2) + 3*x(3))*(1+x(3)^2) + w];
%     f1 = [x(2); -x(1) + x(3); x(1) + (2*x(2) + 3*x(3)) + w];
    loc1 = location(lsupp1, {f1}, p1, 1);

    
        
    %box 2 (control)
    X02 = [];
    
    lsupp2 = loc_support(vars);
    lsupp2 = lsupp2.set_box(L);
    lsupp2.X = [lsupp2.X; r2 >= R0^2];
    lsupp2.X_init = X02;
    lsupp2.Tmax = Tmax;
%     lsupp2.disturb = w^2<=1;
    
    
    % X02 = [t == 0;  (x - [-0.4; 0]).^2 <= 0.01];
    
%     slow2 = 0.75;
    slow2 = 1;
    f2 = [x(2); (-x(1) + x(3)); -x(1) - (2*x(2) + 3*x(3)) + w];
    
    p2 = p1;
%     p2 = [x(1)^2; x(2)^2];

    loc2 = location(lsupp2, f2, p2, 2);
    
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

%BAD BOUNDS
% order = 1; %
% order = 2; %
% order = 3; %
order = 4; %
% order = 5;
%with p1 = -inf
% order = 1%p* = 2.2500000
%     order = 2;  %p* = [0.701374603992632]
%     order = 3;  %p* = [0.473082683617790]
%     order = 4;  %p* = [0.408560944887520]
%     order = 5;  %p* = [0.399423438946940]
%     order = 6;  %p* =[0.389332085206230]


%with p1 = x(1)^2
% order = 1%p* = 2.2500000
%     order = 2;  %p* = 0.701374304679882
%     order = 3;  %p* = 0.472997196248704
%     order = 4;  %p* = 0.408440882926320
%     order = 5;  %p* = 0.400866754178386
%     order = 6;  %p* = 0.390606394334153
    [objective, mom_con, supp_con] =  PM.cons(order);
%     [sol, PM] = PM.run(order, Tmax);    
    sol = PM.run(order) ;
%     fprintf('abs(x1) bound: %0.4f \n', sqrt(sol.obj_rec))
fprintf('bound: %0.4f \n', (sol.obj_rec))
    p_est = sqrt(sol.obj_rec);
    [rr, mm, cc] = PM.recover();
end

%% sample trajectories
if SAMPLE
    rng(25, 'twister')
    smp1 = struct('x', @() sphere_sample(1, 3)'*R0);
%     smp1 = struct('x', @() [R0; 0; 0]);

    smp2 = struct('x', []);
    
    
    LS1 = sampler_uncertain(loc1, smp1);
    
    LS2 = sampler_uncertain(loc2, smp2);
    
    
    LS1.mu = 0.03;
    LS2.mu = 0.03;
    
    HS = sampler_hy({LS1, LS2}, {gfw, gbk});
    
    %     osh = HS.sample_traj(0, [0;0;0.03], 1, 5);
%     Nsample = 1;
%     Nsample = 20;
    Nsample = 50;
%     Nsample = 5;


    [osm, osd] = HS.sample_traj_multi(Nsample, Tmax);
    
    t_end = cellfun(@(o) o.t_end, osm);
    
end

%% plot trajectories
if PLOT
    CPlot = cube_plotter(osm, osd, R0, R1);
    CPlot.cube_plot(sqrt(sol.obj_rec));
    CPlot.nonneg_loc();
    CPlot.aux_plot();
    CPlot.nonneg_jump();

%     CPlot.nonneg_jump();
%     CPlot.objective_plot();
end