classdef pend_manager
    %PEND_MANAGER Formulates and solves the peak estimation problem for
    %stabilization of a pendulum at the angle theta=pi. Also samples
    %trajectories
    %   Detailed explanation goes here
    
    properties
%         param = struct('k', 1, 'delta', 1, 'epsilon', 0.1, 'S'
        
        vars = struct('t', 0, 'x', 0);

        k = 1;          %swing up gain
        delta = 1;      %LQR determination
        epsilon = 0.1;  %energy shaping determination
        eps_tol = 0.01; %wait and LQR are supported at E = +-(epsilon + eps_tol)
        Q = [2 0; 0 1]; %LQR state penalty
        R = 5;          %LQR input penalty
        
        S;              %LQR positive definite matrix
        K;              %LQR input u = K x
         
        PM;             %peak manager
    end
    
    methods
        function obj = pend_manager(k, delta, epsilon, Q, R)
            %PEND_MANAGER Construct an instance of this class
            %   Detailed explanation goes here
            if nargin
                obj.k = k;
                obj.delta = delta;
                obj.epsilon = epsilon;
                obj.Q = Q;
                obj.R = R;
            end
            
            %get the LQR gain and parameters
            A = [0 1; 1 0];
            B = [0; -1];
            [obj.K, obj.S, ~] = lqr(A, B, obj.Q, obj.R);
            
            %make the variables
            mset clear;
            mpol('t', 1, 1);
            mpol('x', 3, 1);
            obj.vars = struct('t', t, 'x', x);
        end        
        
           %% Peak Manager
        function obj = make_manager(obj,  X0, partitions, Tmax, w_max)
            %MAKE_MANAGER form the peak estimation manager for the pendulum
            %system. This includes setting up locations and guards with the
            %height partition 
            %
            %Inputs:
            %X0:    initial support set in terms of variables or a scalar
            %Tmax:  total time of peak analysis
            %w_max: maximum angular velocity of pendulum
            %partitions:    number of height partitions
            %
            %Outputs:
            %PM:    An instance of peak_manager_hy for the pendulum           
            
            if nargin < 4
                Tmax = 10;
            end

            if nargin < 5
                w_max = 2.5;
            end
            
            if nargin < 3
                partitions = 0;
            end
            
            %form the variables

            
            nodes = obj.cheb_nodes(partitions);
            
            loc = cell(partitions+1, 1);
            guard_same_height = cell(partitions+1, 1);
            
            %define locations for each height partition
            id0 = [1;1];
            for i = 1:(partitions+1)
                cmin = nodes(i);
                cmax = nodes(i+1);
                
                [loc_curr, guard_curr] = obj.make_locations(id0, X0, Tmax, w_max, cmin, cmax);
                
                loc{i} = loc_curr;
                id0 = id0 + [4; length(guard_curr)];
                guard_same_height{i} = guard_curr;
            end
            
            %two-way guards between different heights
            guard_diff_height = cell(partitions, 1);
            for i = 1:(partitions)
%                 guard_curr = obj.make_guards_diff_height(id0, loc{i}, loc{i+1}, w_max, nodes(i+1));
                guard_curr = obj.make_guards_diff_height_merged(id0, loc{i}, loc{i+1}, w_max, nodes(i+1));
%                 guard_curr = obj.make_guard_diff_height();
                guard_diff_height{i}= guard_curr;
                id0 = id0 + [0; length(guard_curr)];
            end
            
            %concatenate the locations together 
            loc_cell = cellfun(@(l) {l.up; l.down; l.wait; l.lqr}, loc, 'UniformOutput', false);
            locations = cat(1, loc_cell{:});
            
            %concatenate the guards together 
            guards = [cat(1, guard_same_height{:}); cat(1, guard_diff_height{:})];
            
            obj.PM = peak_manager_hy(locations, guards);
            
            
        end
        
        function [locs, guards] = make_locations(obj, id0, X0, Tmax, w_max, cmin, cmax)
            %MAKE_LOCATIONS form locations for the pendulum system between
            %heights (cosines) [cmin, cmax]. 
            %
            %Inputs:
            %id0:   Index for locations and guards [id_loc, id_guard]
            %X0:    Initial support set
            %Tmax:  total time of peak analysis
            %w_max: maximum angular velocity of pendulum
            %cmin:  minimal cosine (maximum height)
            %cmax:  maximal cosine (minimum height)
            %
            %Output:
            %locs:  A data structure containing the four locations 
            %           (up, down, wait, lqr)
            %guards:a cell array containing the guards between the
            %       locations
            
            %split off variables
            t = obj.vars.t;
            x = obj.vars.x;
            c = obj.vars.x(1);
            s = obj.vars.x(2);
            w = obj.vars.x(3);
            
            %objective to maximize
            objective = x(3)^2;        %maximum angular velocity

            %dynamics
            E = 0.5*w^2 - c;
%             f_wait = Tmax*[-s*w; c*w; -s];
            f_wait = [-s*w; c*w; -s];
            
%             f_lqr = f_wait + Tmax*[0; 0; 1]*obj.u_lqr(x);            
            f_lqr = f_wait + [0; 0; 1]*obj.u_lqr(x);            
            quad = obj.quad_lqr(x);
            
%             f_swing = f_wait + Tmax*[0; 0; 1]*obj.u_swing(x);
            f_swing = f_wait + [0; 0; 1]*obj.u_swing(x);

            
            %support sets
%             c >= cmin; c <= cmax; %may be more numerically stable
            X = [w^2 <= w_max^2;
                 c^2+s^2 == 1; (c - cmin)*(cmax - c) >= 0];
            %location support sets
            X_up = [X; E <= 1-obj.epsilon];   %add energy
            X_down = [X; E >= 1+obj.epsilon]; %reduce energy
            X_wait = [X; E <= 1+obj.epsilon; E >= 1- obj.epsilon; quad >= obj.delta];
            X_lqr = [X; E <= 1+obj.epsilon; E >= 1- obj.epsilon; quad <= obj.delta];

            
            %initial sets            
            if isnumeric(X0)
                %single scalar value as initial support
                %figure out where X0 is located among locations
                X0_up = []; X0_down = []; X0_wait = []; X0_lqr = [];
                X0 = [0; obj.trig_lift(X0)];
                if eval(X_up, [t; x], X0)
                    X0_up = X0; %x0 is in swing-up location
                elseif eval(X_down, [t; x], X0)
                    X0_down = X0; %x0 is in swing-down location
                elseif eval(X_wait, [t; x], X0)
                    X0_wait = X0; %x0 is in waiting location
                elseif eval(X_lqr, [t; x], X0)
                    X0_lqr = X0; %x0 is in lqr location
                end
            else
                %X0 is a support constraint @supcon
                %not sure about any further pruning at this point. 
                X0 = [t == 0; c^2+s^2 == 1; (c - cmin)*(cmax - c) >= 0; X0];
                X0_up = X0; X0_down = X0; X0_wait = X0; X0_lqr = X0;
            end
            
            %assign locations
            id_loc = id0(1);
            id_g   = id0(2);
            
            param = struct('c_range', [cmin, cmax], 'epsilon', obj.epsilon,...
                           'delta', obj.delta, 'S', obj.S, 'w_max', w_max, ...
                           'K', obj.K, 'k', obj.k, 'eps_tol', obj.eps_tol, 'Tmax', Tmax);
            
%             loc_up   = pend_location(id_loc+0, obj.vars, X_up, X0_up, f_swing, param, pend_type.UP);
            loc_up   = pend_location_up(id_loc+0, obj.vars, X_up, X0_up, f_swing, param);
            loc_down = pend_location_down(id_loc+1, obj.vars, X_down, X0_down, f_swing, param);
            loc_wait = pend_location_wait(id_loc+2, obj.vars, X_wait, X0_wait, f_wait, param);
            loc_lqr  = pend_location_lqr(id_loc+3, obj.vars, X_lqr, X0_lqr, f_lqr, param);
            
            %guards of current height range
            %bug over here, wrong signs on energy criterion
            X_up_wait = [X; E == 1 - obj.epsilon; quad >= obj.delta];
            X_up_lqr = [X; E == 1 - obj.epsilon; quad <= obj.delta];
            X_down_wait = [X; E == 1 + obj.epsilon; quad >= obj.delta];
            X_down_lqr = [X; E == 1 + obj.epsilon; quad <= obj.delta];
            X_wait_lqr = [X; E <= 1 + obj.epsilon; E >= 1 - obj.epsilon; quad == obj.delta];

            Reset = x;%identity reset map
            g_up_wait = guard(id_g+0, obj.vars, loc_up, loc_wait, X_up_wait, Reset);      
            g_up_lqr = guard(id_g+1, obj.vars, loc_up, loc_lqr, X_up_lqr, Reset);
            g_down_wait = guard(id_g+2, obj.vars, loc_down, loc_wait, X_down_wait, Reset);            
            g_down_lqr = guard(id_g+3, obj.vars, loc_down, loc_lqr, X_down_lqr, Reset);
            g_wait_lqr = guard(id_g+4, obj.vars, loc_wait, loc_lqr, X_wait_lqr, Reset);


            locs = struct('up', loc_up, 'down', loc_down, ...
                'wait', loc_wait, 'lqr', loc_lqr);
            
            loc_up.g.wait = g_up_wait; loc_up.g.lqr = g_up_lqr;
            loc_down.g.wait = g_down_wait; loc_down.g.lqr = g_down_lqr;
            loc_wait.g.lqr = g_wait_lqr; 
            
            guards = {g_up_wait; g_up_lqr; g_down_wait; g_down_lqr; g_wait_lqr};
        end
                     
        
        function guards = make_guards_diff_height_merged(obj, id0, loc1, loc2, w_max, c_join)            
            %MAKE_GUARDS_DIFF_HEIGHT_MERGED form guard measures between 
            %locations with adjacent heights. Still the identity reset map
            %MERGED: measures with s = +-sin(c_join).
            %
            %Inputs:
            %id0:   Index for locations and guards [id_loc, id_guard]
            %loc1, loc2: structs containing adjacent locations
            %Tmax:  total time of peak analysis
            %w_max: maximum angular velocity of pendulum
            %c_join:common height between location (cosine_join)            
            %
            %Output:
            %guards:a cell array containing the guards between the
            %       locations joining differnet heights
            
            %identity reset map
            
            %split off variables
            t = obj.vars.t;
            x = obj.vars.x;
            c = obj.vars.x(1);
            s = obj.vars.x(2);
            w = obj.vars.x(3);          
            
            %find sin coordinates
%             s_join = sqrt(1-c_join^2);
            
            %evaluate separators
            %q = [-s; w]'*obj.S*[-s; w];
            Reset = x;
            E = 0.5*x(3)^2 - c_join;                
%             quad = obj.quad_lqr(x);
            X = [t*(1-t) <= 1; w^2 <= w_max^2; c==c_join; s^2 == 1 - c_join^2];
            
            %positive root of sin

            quad =  [-s; w]'*obj.S*[-s; w];
            

            %location support sets
            X_up = [X; E <= 1-obj.epsilon];   %add energy
            X_down = [X; E >= 1+obj.epsilon]; %reduce energy
            X_wait = [X; E <= 1+obj.epsilon; E >= 1- obj.epsilon; quad >= obj.delta];
            X_lqr = [X; E <= 1+obj.epsilon; E >= 1- obj.epsilon; quad <= obj.delta];
                                   
            %BAD COPY PASTING HERE form the guards
            %somewhere down the line prune out heights where there is no
            %lqr point, the lqr measure will not be needed
            id_g = id0(2);
            g_up   = obj.height_guard(id_g+0, loc1.up, loc2.up, X_up, Reset);
            g_down = obj.height_guard(id_g+2, loc1.down, loc2.down, X_down, Reset);
            g_wait = obj.height_guard(id_g+4, loc1.wait, loc2.wait, X_wait, Reset);
            g_lqr  = obj.height_guard(id_g+6, loc1.lqr, loc2.lqr, X_lqr, Reset);

            
            %pack up the guards
            guards = [g_up; g_down; g_wait; g_lqr];
        end
        
%         function guards = make_guards_diff_height(obj, id0, loc1, loc2, w_max, c_join)            
%             %MAKE_GUARDS_DIFF_HEIGHT form guard measures between locations
%             %with adjacent heights. Still the identity reset map
%             %
%             %Inputs:
%             %id0:   Index for locations and guards [id_loc, id_guard]
%             %loc1, loc2: structs containing adjacent locations
%             %Tmax:  total time of peak analysis
%             %w_max: maximum angular velocity of pendulum
%             %c_join:common height between location (cosine_join)            
%             %
%             %Output:
%             %guards:a cell array containing the guards between the
%             %       locations joining differnet heights
%             
%             %identity reset map
%             
%             %split off variables
%             t = obj.vars.t;
%             x = obj.vars.x;
%             c = obj.vars.x(1);
%             s = obj.vars.x(2);
%             w = obj.vars.x(3);          
%             
%             %find sin coordinates
%             s_join = sqrt(1-c_join^2);
%             
%             %evaluate separators
%             %q = [-s; w]'*obj.S*[-s; w];
%             Reset = x;
%             E = 0.5*x(3)^2 - c_join;                
% %             quad = obj.quad_lqr(x);
%             X_base = [t*(1-t) <= 1; w^2 <= w_max^2; c==c_join];
%             
%             %positive root of sin
% 
%             quad_pos =  [-s_join; w]'*obj.S*[-s_join; w];
%             X_pos =  [X_base; s == s_join];
% 
%             %location support sets
%             X_up_pos = [X_pos; E <= 1-obj.epsilon];   %add energy
%             X_down_pos = [X_pos; E >= 1+obj.epsilon]; %reduce energy
%             X_wait_pos = [X_pos; E <= 1+obj.epsilon; E >= 1- obj.epsilon; quad_pos >= obj.delta];
%             X_lqr_pos = [X_pos; E <= 1+obj.epsilon; E >= 1- obj.epsilon; quad_pos <= obj.delta];
%             
%             %negative root of sin
%             quad_neg =  [s_join; w]'*obj.S*[s_join; w];            
%             X_neg =  [X_base; s == -s_join];
%             
%             X_up_neg = [X_neg; E <= 1-obj.epsilon];   %add energy
%             X_down_neg = [X_neg; E >= 1+obj.epsilon]; %reduce energy
%             X_wait_neg = [X_neg; E <= 1+obj.epsilon; E >= 1- obj.epsilon; quad_neg >= obj.delta];
%             X_lqr_neg = [X_neg; E <= 1+obj.epsilon; E >= 1- obj.epsilon; quad_neg <= obj.delta];
%             
%             %BAD COPY PASTING HERE form the guards
%             %somewhere down the line prune out heights where there is no
%             %lqr point, the lqr measure will not be needed
%             id_g = id0(2);
%             g_up_pos   = obj.height_guard(id_g+0, loc1.up, loc2.up, X_up_pos, Reset);
%             g_down_pos = obj.height_guard(id_g+2, loc1.down, loc2.down, X_down_pos, Reset);
%             g_wait_pos = obj.height_guard(id_g+4, loc1.wait, loc2.wait, X_wait_pos, Reset);
%             g_lqr_pos = obj.height_guard(id_g+6, loc1.lqr, loc2.lqr, X_lqr_pos, Reset);
%             
%             g_up_neg   = obj.height_guard(id_g+8, loc1.up, loc2.up, X_up_neg, Reset);
%             g_down_neg = obj.height_guard(id_g+10, loc1.down, loc2.down, X_down_neg, Reset);
%             g_wait_neg = obj.height_guard(id_g+12, loc1.wait, loc2.wait, X_wait_neg, Reset);
%             g_lqr_neg  = obj.height_guard(id_g+14, loc1.lqr, loc2.lqr, X_lqr_neg, Reset);
%             
%             %pack up the guards
%             guards = [g_up_pos; g_down_pos; g_wait_pos; g_lqr_pos;
%                       g_up_neg; g_down_neg; g_wait_neg; g_lqr_neg];
%         end
        
        function g_out = height_guard(obj, id, loc_fw, loc_bk, X, Reset)
            g_fw = guard(id, obj.vars, loc_fw, loc_bk, X, Reset);
            g_bk = guard(id+1, obj.vars, loc_bk, loc_fw, X, Reset);
            g_out = {g_fw; g_bk};
            
            loc_fw.g.raise = g_fw;
            loc_bk.g.lower = g_bk;
        end
        
        %% helper functions
        
        function x_lift = trig_lift(obj, x_angle)
            %lift from angle to semialgebraic representation
            x_lift = [cos(x_angle(1, :)); sin(x_angle(1, :)); x_angle(2, :)];
        end
        
        function E = energy(obj, x)    
            %energy function in lifted coordinates
            E = 0.5*x(3, :).^2 - x(1, :);
        end
        
        function E = energy_ang(obj, x_angle)
            E = obj.energy(obj.trig_lift(x_angle));
        end
        
        function q = quad_lqr(obj, x)
            %LQR ellipsoid in lifted coordinate (cylinder in s-w)
            s = x(2);
            w = x(3);
            q = [-s; w]'*obj.S*[-s; w];
        end
                
        function q = quad_lqr_ang(obj, x_angle)
            q = obj.quad_lqr(obj.trig_lift(x_angle));
        end
                
        
        function [obj, w_max_est, sol] = run(obj, order)
            %SOLVE solve the maximum angular velocity estimation problem
            if nargin < 2
                order = 2;
            end
            
            % [objective_mom, mom_con, supp_con] =  PM.peak_cons(order);
            [obj.PM, sol] = obj.PM.run(order);
%             [sol, obj] = run@peak(order);
            w_max_est = sqrt(sol.obj_rec);
            
            
        end
        
        %% Controllers
        function u = u_lqr(obj, x)
            u = obj.K*[-x(2); x(3)];
%             K*[-sin(x(1)); x(2)];
        end
        
        function u = u_swing(obj, x)
            u = -obj.k * x(3) * (obj.energy(x)-1);
        end
        
        function u = u_lift(obj, x)
            E = obj.energy(x);
            if (E <= (1-obj.epsilon)) || (E >= (1 + obj.epsilon))
%                 u = obj.u_swing(x);
                u = -obj.k * x(3) * (E-1);
            else
                q = obj.quad_lqr(x);
                if q <= obj.delta
%                     u = obj.u_lqr(x);
                    u = obj.K*[-x(2); x(3)];
                else
                    u = 0;
                end
            end
        end
        
        function u = u_ang(obj, x)
            %hybrid controller
            x_lift = obj.trig_lift(x);
            E = obj.energy(x_lift);
            if (E <= (1-obj.epsilon)) || (E >= (1 + obj.epsilon))
                u = obj.u_swing(x_lift);
            else
                q = obj.quad_lqr(x_lift);
                if q <= obj.delta
                    u = obj.u_lqr(x_lift);
                else
                    u = 0;
                end
            end
        end                
        
        function xdot = f_ang(obj, t, x)
            %main dynamics
%             u = obj.u_ang(x);
            x_lift = obj.trig_lift(x);
            u = obj.u_lift(x_lift);
            xdot = [x_lift(3); -x_lift(2) +  u];
        end
       
        
        %% Partitions        
        function c_list = cheb_nodes(obj, N)
            %return partition of [-1, 1] by Chebyshev nodes
            r = 2*(1:N);
            c_list = [-1 -cos(pi/(2*N) * (r - 1)) 1];
        end
        
        function c_list = std_nodes(obj, N)
            %return partition of [-1, 1] by linear spacing
            c_list = linspace(-1, 1, N+2);
%             c_list = cos(pi/(2*N) * (r - 1));
        end
        
        function [x_lqr, theta_lqr] = max_lqr(obj, ang_max)
            %MAX_LQR Find the point (angle) that is at the extreme range of
            %LQR control. After that angle, the only locations are
            %swing-up, swing-down, or wait.
            
            if nargin < 2
                ang_max = 5;
            end
            
            xs = sdpvar(2, 1); %[sin(theta), w]
            xs_lim = [1; ang_max];
            L = chol(obj.S);
            F = [norm(L*[-xs(1); xs(2)], 2) <= sqrt(obj.delta);
                 xs <=  xs_lim; 
                 xs >= -xs_lim];
            h = -xs(1);

            sol = optimize(F, h, sdpsettings('solver', 'mosek', 'verbose', 0));

            % double(xs)

            sx = double(xs(1));
            cx = -sqrt(1-sx^2);
            theta_lqr = atan2(sx, cx)*180/pi;
            
            x_lqr = [cx; double(xs)];
        end
        
     
        
        %% Simulation
        function out_sim = pend_sim(obj, x0, Tmax)
            %based on peak/utils/pend_sim
            %assume time starts at 0
            options = odeset('RelTol',1e-9, 'MaxStep', 0.01);
            [time_accum, x_accum] = ode15s(@(t, x) obj.f_ang(t, x), [0, Tmax], x0, options);
            
            out_sim = struct;
            out_sim.t = time_accum;
            out_sim.theta = x_accum;
            out_sim.Tmax = Tmax;
            
            x_trig = [cos(x_accum(:, 1)), sin(x_accum(:, 1)), x_accum(:, 2)];
            out_sim.x = x_trig;
            
            %next step: hybrid trajectories based on locations
        end
        
        
%         function out_sim = pend_sim_loc(obj, t0, x0, Tmax)
%             %PEND_SIM_LOC simulate the pendulum at the current location id0
%         end
        
        function [out_sim] = pend_sim_nonneg(obj, x0, Tmax)
            %based on peak_manager_hy/sample_traj
            %sample a trajectory and evaluate nonnegative functions
            %necessarily more complicated than obj.pend_sim()
            t_curr = 0;
            x_curr = x0;
            xl_curr = obj.trig_lift(x_curr);
                
            id_curr = find(obj.PM.supp_loc_eval(0, x0));
            
            out_sim.sim = {};
            out_sim.jump = {};
            zeno_count = zeros(length(obj.PM.guards), 1);
            
            %simulate the trajectory

            
            while (t_curr < Tmax)
                %sample within location
                loc_curr = obj.PM.loc{id_curr};
%                 event_curr = @(t, x) obj.PM.loc_event(t/Tmax, x, id_curr); 
                event_curr = @(t, x) loc_curr.supp_event(t, x);
                curr_ode_options = odeset('Events',event_curr, 'RelTol', 1e-7, 'AbsTol', 1e-8,...
                                      'MaxStep', 0.01);
                out_sim_curr = struct;
                
                %dynamics in each location, should be faster
                f_curr = @(t, x) loc_curr.f_eval(t, x);
                [time_accum, x_accum] = ode15s(f_curr, [t_curr, Tmax], x_curr, curr_ode_options);
                
                %store trajectory in location                                                               
                out_sim_curr.t = time_accum;
                out_sim_curr.x = x_accum;
                out_sim_curr.xl = obj.trig_lift(x_accum')';
                out_sim_curr.E = obj.energy_ang(x_accum')';
                t_curr = out_sim_curr.t(end);
                x_curr = out_sim_curr.x(end, :)';                
                xl_curr = obj.trig_lift(x_curr)';
                
                if loc_curr.dual.solved
                    out_sim_curr.nonneg = loc_curr.nonneg_eval(time_accum'/Tmax, out_sim_curr.xl')';
                    out_sim_curr.v = loc_curr.v_eval(time_accum'/Tmax, out_sim_curr.xl')';
                end
                out_sim_curr.objective = x_accum(:, 2).^2;
                out_sim_curr.id = id_curr;
                out_sim_curr.type = loc_curr.type;
                out_sim_curr.u = loc_curr.u_eval_trig(out_sim_curr.xl);
                
                out_sim.sim{end+1} = out_sim_curr;
                


                %figure out jump
                %a hack to use x_prev since the event location is
                %incompatible 
                poss_g = loc_curr.choose_guard(t_curr, x_curr);
                if  ~isempty(poss_g)
                    g_new = poss_g;
                    zeno_count(g_new.id) = zeno_count(g_new.id) + 1;
                    Rx = x_curr;
                    %identity reset map
                    
                    %track the nonnegativity  
                    jump_curr = struct('t', t_curr, 'x', x_curr, 'x_jump', ...
                        Rx, 'guard', g_new.id);
                    
                    if g_new.dual.solved
                        %TODO
                        jump_curr.nonneg = g_new.nonneg(t_curr/Tmax, obj.trig_lift(x_curr));
                    end
                    out_sim.jump{end+1} = jump_curr;
                    
                    
                    %complete the jump
                    id_curr = g_new.dest.id;
                    x_curr = Rx;
                    
                    if zeno_count(g_new.id) > g_new.zeno_cap
                        %maximum number of jumps is exceeded
                        break
                    end
                else
                    %no available guards to jump to
                    %outside support, end of trajectory
                    break
                end
            end
            out_sim.t_end = t_curr;
        end
        
        
        function [out_sim_multi] = pend_sample_multi(obj, init_sampler, Tmax)
            if isnumeric(init_sampler.init)
                %given sample points
                N = size(init_sampler.init, 2);
                out_sim_multi = cell(N, 1);               
                %parallel code requires splitting off separate objects
                for i = 1:N                    
                    x0 = init_sampler.init(:, 2);                    
                    out_sim_multi{i} = obj.pend_sim_nonneg(x0, Tmax);
                end
                
            else
                %random sample.
                N = init_sampler.N;
                out_sim_multi = cell(N, 1);
                for i = 1:N                    
                    [x0] = init_sampler.init();                                        
                    out_sim_multi{i} = obj.pend_sim_nonneg(x0, Tmax);
                end
            end                                                
            
        end
        
        function out_sim_multi = nonneg_eval(obj, out_sim_multi, Tmax)
        %evaluate nonnegative trajectories again after a new run
        
            for i = 1:length(out_sim_multi)
                %nonneg of locations
                for j = 1:length(out_sim_multi{i}.sim)
                    loc_curr = obj.PM.loc{out_sim_multi{i}.sim{j}.id};
                    out_sim_curr = out_sim_multi{i}.sim{j};
                    t_curr = out_sim_curr.t;
                    xl_curr = out_sim_curr.xl;
                    if loc_curr.dual.solved
                        out_sim_multi{i}.sim{j}.nonneg = loc_curr.nonneg_eval(t_curr'/Tmax, xl_curr')';
                        out_sim_multi{i}.sim{j}.v = loc_curr.v_eval(t_curr'/Tmax, xl_curr')';
                    end
                end
                
                %nonneg of guards
                for j = 1:length(out_sim_multi{i}.jump)
                    g_curr = obj.PM.guards{out_sim_multi{i}.jump{j}.guard};
                     if g_curr.dual.solved
                        %TODO
                        t_curr = out_sim_multi{i}.jump{j}.t;
                        xl_curr = obj.trig_lift(out_sim_multi{i}.jump{j}.x);
                        out_sim_multi{i}.jump{j}.nonneg = g_curr.nonneg(t_curr/Tmax, xl_curr);
                    end
                end
            end
            
        end
        
        
        %% initial samplers
        
        function [x_sample, xl_sample] = pend_sample_point_box(obj, c_range, w_range)
           c_sample = c_range(1) + rand()*diff(c_range); 
           w_sample = w_range(1) + rand()*diff(w_range);
           s_mag = 1 - c_sample^2;
           s_sample = sign(2*rand()-1)*s_mag;
           
           x_sample = [atan2(s_sample, c_sample); w_sample];
           xl_sample = [c_sample; s_sample; w_sample];
           
        end
        
        function [x_sample, xl_sample] = pend_sample_circ(obj, w_range)
            w_sample = w_range(1) + rand()*diff(w_range);
            th = 2*pi*rand();
            x_sample = [th; w_sample];
            xl_sample = [cos(th); sin(th); w_sample];
        end
        
        
        function [x_sample, xl_sample] = pend_sample_point_ball(obj, c_range, w_range)

           ball_raw = ball_sample(1, 2);

           c_sample = c_range(1) + ball_raw(1)*diff(c_range); 
           w_sample = w_range(1) + ball_raw(1)*diff(w_range);
           s_mag = 1 - c_sample^2;
           s_sample = sign(2*rand()-1)*s_mag;
           
           x_sample = [atan2(s_sample, c_sample); w_sample];
           xl_sample = [c_sample; s_sample; w_sample];
           
        end
        

    end
end

