function [event_eval, terminal, direction] = stoch_event_fw(t, x, param)
    %event function for @ode15 or other solver
    %Start control of the system (it is too far away from the origin)
    Npt = size(x, 2);
    event_eval = zeros(1, Npt);
    for i = 1:Npt
        xcurr = x(:, i);
        tcurr = t(:, i);               
        
%         crit_bk = sum(xcurr.^2);
        crit_fw = param.K*x(1)^2 + x(2)^2 + x(3)^2;

        event_eval(i) = 2*all([crit_fw <= param.R1])-1; 
    end

    %stop integrating when the system falls outside support

    terminal = 1;
    direction = 0;                        
end