function [osd] = rl_sampler(x0, Tmax, zeno_count)
%RL_SAMPLER Sample the right-left-wrap routine
%I am having a terrible amount of numerical difficulties with tolerances in
%sampler_hy regarding guards. I am therefore using a dedicated sampler
%here.

R = 1;
f = @(x) [-x(2) + x(1)*x(2) + 0.5;
         -x(2) - x(1) + x(1)^3];

event = @(x) all(abs(x) <= R);

osd = struct('locations', {}, 'guards', {{}, {}});

tcurr = 0;

end

