classdef rl_location< location
    %RL_LOCATION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Access=private)
        R = 1;
    end 
    
    methods
        function obj = rl_location(lsupp, f, p, loc_id)
            %BALL2D_FLAT_LOC Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 4
                loc_id = 1;
            end
            obj@location(lsupp, f, p, loc_id);            
        end
        
        function supp_out = supp_eval(obj, t ,x_in)
            %support evaluation of the bouncing ball in a curved basin
            %is the bottom of the basin, so y <= 0
            %will untangle later
            supp_out = all(abs(x_in) <= obj.R);
        end

        function x_proj = supp_proj(obj, x_in)
            %project onto the support set (box)
            mask = abs(x_in) > obj.R;

            x_proj = x_in;

            x_proj(mask) = obj.R*sign(x_in(mask));
        end
        
    end
end

