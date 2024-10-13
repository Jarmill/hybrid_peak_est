%newest script to formulate pendulum swingup experiment
SETUP = 1;
SOLVE = 1;
SAMPLE = 0;
EVAL = 1;
PLOT = 1;
if SETUP
%there was an optimal setting. What was it?
%maybe that was with starting at pi
    

%set partitions and orders
Npartitions = 0;
% order = 1; %p* = 4.2386 (valid)
% order = 2; %p* = 4.2006 (valid)
% order = 3; %p* = 4.2000 (valid)
order = 4; %p* = 4.1998
% order = 5; %p*=4.1885 (invalid)
% order = 6; %p*= 3.6948 (invalid)

% Npartitions = 3;
% order = 2; %p* = 4.2005
% order = 3; %p* = 4.1997




%the experiment

Tmax = 10;
    


w_range_0 = 0.5*[-1; 1];

    PM = pend_manager();
%     supp0 = x0;
%     supp0 = [PM.vars.t == 0; PM.vars.x == PM.trig_lift(x0)];
    supp0 = [PM.vars.x(3)^2 <= w_range_0(2)^2];
    PM = PM.make_manager(supp0, Npartitions, Tmax);
end
if SOLVE

orderlist = 1:6;
    p_order = zeros(length(orderlist), 1);
    time_order = zeros(length(orderlist), 1);
     for i = 1:length(orderlist)
    
    %     [objective, mom_con, supp_con, len_dual] =  PM.cons(order);
      % tic
            [PM, w_max_est, sol] = PM.run(orderlist(i));
            time_order(i) = sol.solver_time;
            p_order(i) = sol.obj_rec;
    %     [sol, = PM.run(order) ;
    %     fprintf('abs(x1) bound: %0.4f \n', sqrt(sol.obj_rec))
    fprintf('maximum w(3)^2: %0.4f', sol.obj_rec);
        p_est = sqrt(sol.obj_rec);
        % [rr, mm, cc] = PM.recover();
    
        obj_rec = sol.obj_rec;
        save('pendulum_experiment.mat', 'time_order', 'p_order', 'orderlist');
     end
% 
%     [PM, w_max_est, sol] = PM.run(order);
% %     [sol, PM] = PM.run(order);
%     [op, mom_out, corner_out] = PM.PM.recover();
%     op_all = all(op);
    
    %without partitions at order 2:
    %w_max_est = 2.049354326410954

    %2 partitions at order 2
    %w_max_est = 1.901106516822999
end
%% sample
if SAMPLE
        
%     osn = PM.pend_sim_nonneg(x0, 10);


    rng(10, 'twister')
    Nsample = 250;
    ps_func = @() PM.pend_sample_point_box([-1; 1], w_range_0);
    ps = struct('N', Nsample, 'init', ps_func);
    osm = PM.pend_sample_multi(ps, Tmax);

end

%% evaluate nonneg
if EVAL
    osm2 = PM.nonneg_eval(osm, Tmax);
end

%% plot
if PLOT

    PP = pend_plotter(osm2);
    PP.state_plot_trig(w_range_0, w_max_est);
    PP.nonneg_jump();
    PP.nonneg_loc();
    PP.nonneg_loc();
    PP.cost_plot(w_max_est);
    PP.track_plot_trig();    
    PP.aux_plot();
end 
