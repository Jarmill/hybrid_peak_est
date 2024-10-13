classdef pend_location_lqr < pend_location
    %PEND_LOCATION_LQR Pendulum location in the LQR phase
    %quasi-linear control to stabilize the pendulum
    
    properties(Access=private)
        epsilon_;
        w_max_;
        c_range_;
        delta_;
        S_;
        K_;
    end
    
    
    methods
        function obj = pend_location_lqr(id, vars, supp, supp0, f, param)
            %PEND_LOCATION_LQR Construct an instance of this class
            %   Detailed explanation goes here            
            obj@pend_location(id, vars, supp, supp0, f, param, pend_type.LQR);
            obj.epsilon_ = obj.param.epsilon + obj.param.eps_tol;
            obj.w_max_ = obj.param.w_max;
            obj.delta_ = obj.param.delta;
            obj.S_ = obj.param.S;
            obj.K_ = obj.param.K;
            obj.c_range_ = obj.param.c_range;
%             obj.eps_tol_ = obj.param.eps_tol;
        end
        
        function supp_out = supp_eval(obj, t ,x)
            %support evaluation in the LQR location
            
            c = cos(x(1));
            s = sin(x(1));
            w = x(2);
            w_out = (w <= obj.w_max_) && (w >= -obj.w_max_);
            
            c_out = (c >= obj.c_range_(1)) && (c <= obj.c_range_(2));
            E = 0.5*w^2 - c;
            E_out = (E >= (1 - obj.epsilon_)) && (E <= (1 + obj.epsilon_));
            
            quad = [-s w]*obj.S_*[-s; w];
            q_out = (quad <= obj.delta_);
            
            supp_out = c_out && E_out && q_out && w_out;                        
        end    
        
        function u = u_eval(obj, x)                 
            %evaluate the control action
            s = sin(x(:, 1));
            w = x(:, 2);
            u = [-s, w]*obj.K_';
        end
        
        function u = u_eval_trig(obj, x)                 
            %evaluate the control action
            s = x(:, 2);
            w = x(:, 3);
            u = [-s, w]*obj.K_';
        end
        
        
        function xdot = f_eval(obj, t, x)
            %evaluate the swing-up dynamics at WAIT
            c = cos(x(1));
            s = sin(x(1));
            w = x(2);   
            
            u = obj.K_*[-s; w];
            
            xdot = [w; -s + u];
        end
        
    end
end

