classdef single_side_sampler < sampler_hy
    %SINGLE_SIDE_SAMPLER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        R = 1; %sidewall
    end
    
    methods
        function obj = single_side_sampler(sampler_in, guards_in, R)
            %SINGLE_SIDE_SAMPLER Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

