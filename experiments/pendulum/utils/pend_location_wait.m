classdef pend_location_wait < pend_location
    %PEND_LOCATION_WAIT Pendulum location in the swing-up phase
    %removing energy from the system    
    
    properties(Access=private)
        epsilon_;
        w_max_;
        c_range_;
        delta_;
        S_;
    end
    
    methods
        function obj = pend_location_wait(id, vars, supp, supp0, f, param)
            %PEND_LOCATION_WAIT Construct an instance of this class
            %   Detailed explanation goes here            
            obj@pend_location(id, vars, supp, supp0, f, param, pend_type.WAIT);
            obj.epsilon_ = obj.param.epsilon + obj.param.eps_tol;
            obj.w_max_ = obj.param.w_max;
            obj.delta_ = obj.param.delta;
            obj.S_ = obj.param.S;
            obj.c_range_ = obj.param.c_range;
            
        end
        
        function supp_out = supp_eval(obj, t ,x)
            %support evaluation in the DOWN location
            
            c = cos(x(1));
            s = sin(x(1));
            w = x(2);
            tol = 1e-2; %a tolerance to prevent chattering in the waiting phase
            w_out = (w <= obj.w_max_) && (w >= -obj.w_max_);
            
            c_out = (c >= obj.c_range_(1)) && (c <= obj.c_range_(2));
            E = 0.5*w^2 - c;
            E_out = ((E+tol) >= (1 - obj.epsilon_)) && ((E-tol) <= (1 + obj.epsilon_));
            
            quad = [-s w]*obj.S_*[-s; w];
            q_out = (quad >= obj.delta_);
            
            supp_out = c_out && E_out && q_out && w_out;
        end     
        
        function u = u_eval(obj, x)                 
            %evaluate the control action
%             s = sin(x(:, 1));
%             w = x(:, 2);
            u = zeros(size(x(:, 1)));
        end
        
        function u = u_eval_trig(obj, x)                 
            %evaluate the control action
%             u = [-s, w]*obj.K_';
            u =  zeros(size(x(:, 1)));
        end
        
        function xdot = f_eval(obj, t, x)
            %evaluate the swing-up dynamics at WAIT
%             c = cos(x(1));
            s = sin(x(1));
            w = x(2);                        
            
            xdot = [w; -s];
        end
        
    end
end

