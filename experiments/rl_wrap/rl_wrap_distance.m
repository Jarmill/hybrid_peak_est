% A single-location case where the left side jumps to the bottom, and the
% right side jumps to the top. Both the top and bottom are purely
% inward-facing.
%
% There is no noise in this example. This script investigates whether the
% auxiliary function v(t, x) decreases along noise-free trajectories



SETUP = 1;
SOLVE = 0;
SAMPLE = 0;
PLOT = 1;
% rl_analyze_traj;


%% setup the locations and guards
if SETUP
    
    mset clear
    mpol('t');
    mpol('x', 2, 1);
    mpol('y', 2, 1);
%     mpol('w', 1, 1);
    w=0;
    vars = struct('t', t, 'x', x, 'y', y);

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
        
    lsupp1 = unsafe_support(vars);
    lsupp1 = lsupp1.set_box(L);
%     lsupp1.X = [lsupp1.X];
    lsupp1.X_init = X01;
    lsupp1.Tmax = Tmax;
%     lsupp1.disturb = w^2<=1;
    
%     p1 = -x'*x;

%     Cu = [-0.5; -0.5];
%     Ru = 0.4;

    Cu = [-0.3; -0.3];
    Ru = 0.4;
    
    theta_c = 5*pi/4;
    w_c = [cos(theta_c); sin(theta_c)];
    c2f = w_c(1)*(y(1) - Cu(1)) + w_c(2) * (y(2) - Cu(2)); 
    
    Xu = [Ru^2 - sum((y-Cu).^2); c2f]>=0;
    lsupp1.X_unsafe = Xu;
    lsupp1.dist = (x-y)'*(x-y);
%     p1 = -sum((x-Cu).^2);

%     p1 = -sum(x);
%     p1 = [];

%     f1 = [-x(2) + x(1)*x(2)/2;
%          -x(2) - x(1) + x(1)^3/4];

    %original
    %k = 0.5;
    k = 0.6;
    f1 = [-x(2) + x(1)*x(2) + k;
         -x(2) - x(1) + x(1)^3];

    CSP = 1;
    if CSP
        loc1 = location_distance_csp(lsupp1, {f1}, 1);
    else
        loc1 = location_distance(lsupp1, {f1}, 1);
    end
     
%     loc1 = rl_location(lsupp1, {f1}, p1, 1);

    
        
  %guards

    %reset maps
    Rleft = [1-x(2); x(1)];     %left->bottom
    Rright = [x(2); x(1)]; %right->top
    
    
    %possible bug: time not defined in support in the reset maps
    Xleft  = [x(2)^2 <= R1^2; x(1)==-R1; t*(1-t)>=0];   
    Xright = [x(2)^2 <= R1^2; x(1)==R1; t*(1-t)>=0];   
    
    gleft  = guard(1, vars, loc1, loc1, Xleft, Rleft);
    gright = guard(2, vars, loc1, loc1, Xright, Rright);
    
    
    Zeno_cap = Inf;
    % Zeno_cap = 5;
%     Zeno_cap = 10;
%     Zeno_cap = 15;
    gleft.zeno_cap = Zeno_cap;
    gright.zeno_cap = Zeno_cap;
end

%% solve the system and get peak estimates
if SOLVE
    
    PM =  distance_manager_hy({loc1}, {gleft, gright});
%     PM =  peak_manager_hy({loc1}, {gleft, gright});

%CSP: distances (norm)

%     Cu = [-0.5; -0.5];
%     Ru = 0.4;
%     order = 1; %0
%     order = 2; %0
%     order = 3; % 0.2157 
%     order = 4; % 0.3295 
%     order=5; %0.3691
%     order = 6; % 0.3761 


    Cu = [-0.3; -0.3];
    order=1; %0
%     order=2; %0
%     order=3; %0.1647
%     order=4; %0.2643
%     order=5; %0.2853
%     order=6; %0.2877
%     Ru = 0.4;
%     
%     theta_c = 5*pi/4;
% order=6;    


orderlist = 1:6;
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
        save('rl_wrap_distance_nocap.mat', 'time_order', 'p_order', 'orderlist');
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
    Nsample = 20;
% Nsample = 30;
%     Nsample = 50;
    % Nsample = 120;
    % Nsample = 
%     Nsample = 5;


    [osm, osd] = HS.sample_traj_multi(Nsample, Tmax);
    
    t_end = cellfun(@(o) o.t_end, osm);
    
end

%% plot trajectories
if PLOT
    load('rl_traj_5.mat', 'osm', 'osd')
    clf
   
    mm = -Inf;
    % dd = @(x) aff_half_circ_dist(x', Ru, theta_c, Cu);

    dd = @(x) -norm(x - [0.5, 0.3]');
    for i = 1:length(osd.locations{1})
        traj = osd.locations{1}{i};
        for j = 1:length(traj.x)
            mm = max(mm, dd(traj.x(j, :)));
        end
    end
    

    RPlot = rl_plotter(osm, osd, R0, C0, Cu);
%     hold on
    RPlot.rl_plot();
    title('hi')
    
%     theta_c = 3*pi/2;
    obj_rec = 0.2877^2;
    theta_half_range = linspace(theta_c-pi/2, theta_c + pi/2, 200);
    Rot_mat = [cos(theta_c+pi/2), -sin(theta_c+pi/2); sin(theta_c+pi/2), cos(theta_c+pi/2)];
    x_dist_eps = Rot_mat*dist_contour(200, Ru, sqrt(obj_rec)) + Cu;
    circ_half = [cos(theta_half_range); sin(theta_half_range)];
    Xu = Cu + circ_half* Ru;
    patch(Xu(1, :), Xu(2, :), 'r', 'Linewidth', 3, 'EdgeColor', 'none', 'DisplayName', 'Unsafe Set')

    plot(x_dist_eps(1, :), x_dist_eps(2, :), 'r', 'DisplayName', 'Distance Contour', 'LineWidth', 2)

%     RPlot.rl_plot(obj_rec);
    
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


function x_dist = dist_contour(Ntheta, R, c)
    %compute a contour at distance c away from the half-circle with N_theta
    %sample points


    theta_q1 = linspace(0, pi/2, Ntheta);
    theta_q2 = linspace(pi/2, pi, Ntheta);
    theta_q34 = linspace(pi, 2*pi, 2*Ntheta);

    %contour level
    

    x_right = [c*cos(theta_q1)+R; c*sin(theta_q1)];
    x_left = [c*cos(theta_q2)-R; c*sin(theta_q2)];
    x_bottom = [(c+R)*cos(theta_q34); (c+R)*sin(theta_q34)];
    x_dist = [x_right, x_left, x_bottom];
end

function dist_out = half_circ_dist(x_in, R)
    %return the L2 distance  between the point x_in and the half circle
    %||x_in||^2 <= R^2 intersect x_in(2) <= 0.
%     reshape(x_in, [], 1);
    if x_in(2) >= 0
        %flat region
        if x_in(1) < -R
            dist_out = hypot(x_in(1)+R, x_in(2));
        elseif x_in(1) > R
            dist_out = hypot(x_in(1)-R, x_in(2));
        else
            dist_out = x_in(2);
        end
    else
        %circle region
        dist_out = max(norm(x_in, 2)-R, 0);
    end

end

function dist_out = aff_half_circ_dist(x_in, R, theta_c, Cu)

    theta_cf = theta_c - 3*pi/2;
    % theta_cf = theta_c;
    Rot_mat = [cos(theta_cf) -sin(theta_cf); sin(theta_cf) cos(theta_cf)];
    x_aff = Rot_mat'*(x_in - Cu);
    
    dist_out = half_circ_dist(x_aff, R);

end