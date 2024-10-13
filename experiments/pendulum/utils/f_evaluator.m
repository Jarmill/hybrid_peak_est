classdef f_evaluator
    %F_EVALUATOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Access=private)
        S;
        K;
        epsilon;
        delta;
    end
    
    methods
        function obj = f_evaluator(S,K, epsilon, delta)
            %F_EVALUATOR Construct an instance of this class
            %   Detailed explanation goes here
            obj.S = S;
            obj.K = K;
            obj.epsilon = epsilon;
            obj.delta = delta;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

