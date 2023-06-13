function [event_eval, terminal, direction] = stoch_event_bk(t, x, param)
    %event function for @ode15 or other solver
    %Stop control of the system (let drift because it is close to the
    %origin)
    Npt = size(x, 2);
    event_eval = zeros(1, Npt);
    for i = 1:Npt
        xcurr = x(:, i);
        tcurr = t(:, i);               
        
        crit_bk = norm(xcurr);

        event_eval(i) = 2*all([abs(x) <= param.L; crit_bk > param.R0])-1; 
    end

    %stop integrating when the system falls outside support

    terminal = 1;
    direction = 0;                        
end