classdef pend_location < location
    %PEND_LOCATION A location for peak estimation of the swing-up pendulum
    %system. Based on the standard utils/location class
    
    properties
        param = struct;
        g = struct('raise', [], 'lower', [], 'wait', [], 'lqr', []);
%         h_idx = 0;  %height index in partition
        type; %type of location (of pend_type)
    end
    
    methods        
        %param = {w_max, Tmax, epsilon, delta, c_range}
        function obj = pend_location(id, vars, supp, supp0, f, param, p_type)
            objective = vars.x(3)^2;
            lsupp = loc_support(vars);
            lsupp.X = supp;
            lsupp.X_init = supp0;
            lsupp.Tmax = param.Tmax;
            obj@location(lsupp, f, objective, id);
            obj.param = param;
            obj.type = p_type;
        end
        
        function supp_out = supp_eval(obj, t ,x)
            %support evaluation of the pendulum system
            c = cos(x(1));
            s = sin(x(1));
            w = x(2);
            w_out = (w <= obj.param.w_max) && (w >= -obj.param.w_max);
            E = 0.5*w^2 - c;
            c_out = (c >= obj.param.c_range(1)) && (c <= obj.param.c_range(2));
            if obj.type == pend_type.UP
                E_out = (E <= (1 - obj.param.epsilon));
                q_out = 1;
            elseif obj.type == pend_type.DOWN
                E_out = (E >= (1 + obj.param.epsilon));
                q_out = 1;
            else
                E_out = (E >= (1 - obj.param.epsilon)) && (E <= (1 + obj.param.epsilon));
                quad = [-s w]*obj.param.S*[-s; w];
                if obj.type == pend_type.LQR
                    q_out = (quad <= obj.param.delta);
                else
                    q_out = (quad >= obj.param.delta);
                end
            end
            
            supp_out = c_out && E_out && q_out && w_out;
                    
        end    
        
        function g_out = choose_guard(obj, t, x)
            %what guard to move to?
            g_out = [];
            c = cos(x(1));
            s = sin(x(1));
            w = x(2);
            
            tol = 1e-8;
            
            if abs(c - obj.param.c_range(1)) <= tol
                %lower the height
                g_out = obj.g.lower;
            elseif abs(c - obj.param.c_range(2)) <= tol
                %raise the height
                g_out = obj.g.raise;
            else
                %transition based on quad and energy
                E = 0.5*w^2 - c;
                E_target = ((E <= 1+obj.param.epsilon) || (E >= 1-obj.param.epsilon));
                quad = [-s w]*obj.param.S*[-s; w];
                if ((obj.type == pend_type.UP) || (obj.type == pend_type.DOWN)) && E_target
                    if quad <= obj.param.delta
                        g_out = obj.g.lqr;
                    else
                        g_out = obj.g.wait;
                    end
                elseif (obj.type == pend_type.WAIT) && (quad <= obj.param.delta)
                    g_out = obj.g.lqr;
                end                                    
            end
        end
       
        
        function x_lift = trig_lift(obj, x_angle)
            %lift from angle to semialgebraic representation
            x_lift = [cos(x_angle(1, :)); sin(x_angle(1, :)); x_angle(2, :)];
        end               
        
    end
end

