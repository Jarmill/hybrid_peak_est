classdef pend_location_up < pend_location
    %PEND_LOCATION_UP Pendulum location in the swing-up phase
    %pumping energy into the system    
    
    properties(Access=private)
        epsilon_;
        w_max_;
        c_range_;
        k_;
    end
    
    methods
        function obj = pend_location_up(id, vars, supp, supp0, f, param)
            %PEND_LOCATION_UP Construct an instance of this class
            %   Detailed explanation goes here            
            obj@pend_location(id, vars, supp, supp0, f, param, pend_type.UP);
            obj.epsilon_ = obj.param.epsilon;
            obj.w_max_ = obj.param.w_max;
            obj.k_ = obj.param.k;
            obj.c_range_ = obj.param.c_range;
        end
        
        function supp_out = supp_eval(obj, t ,x)
            %support evaluation in the UP location
            
            c = cos(x(1));
            s = sin(x(1));
            w = x(2);
            w_out = (w <= obj.w_max_) && (w >= -obj.w_max_);
            
            c_out = (c >= obj.c_range_(1)) && (c <= obj.c_range_(2));
            E = 0.5*w^2 - c;
            E_out = (E <= (1 - obj.epsilon_));
            
            
            supp_out = c_out && E_out && w_out;
        end     
        
        function u = u_eval(obj, x)                 
            %evaluate the control action
            c = cos(x(:, 1));
            w = x(:, 2);
            E = 0.5*w.^2 - c;
%             u = [-s, w]*obj.K_';
            u =  -obj.k_ .* w .* (E-1);
        end
        
        function u = u_eval_trig(obj, x)                 
            %evaluate the control action
            c = x(:, 1);
            w = x(:, 3);
            E = 0.5*w.^2 - c;
%             u = [-s, w]*obj.K_';
            u =  -obj.k_ .* w .* (E-1);
        end
        
        function xdot = f_eval(obj, t, x)
            %evaluate the swing-up dynamics at UP
            c = cos(x(1));
            s = sin(x(1));
            w = x(2);
            
            E = 0.5*w^2 - c;
            
            u =  -obj.k_ * w * (E-1);
            
            xdot = [w; -s + u];
        end
        
    end
end

