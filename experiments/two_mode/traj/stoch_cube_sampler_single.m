function [osd] = stoch_cube_sampler_single(x0, Tmax, param)
%stoch_cube_sampler_single Sample the stochastic cube routine
%I am having a terrible amount of numerical difficulties with tolerances in
%sampler_hy regarding guards. I am therefore using a dedicated sampler
%here.

%based on sampler_hy.sample_traj

%param: stuct with fields
%   R0:     inner radius
%   R1:     outer radius
%   K:      ellipicity of outer set
%   sigma:  influence of Brownian noise
%   dt:     time spacing (1e-3)

%dynamics
f1 = @(t, x) [x(2); (-x(1) + x(3)); x(1) + (2*x(2) + 3*x(3))*(1+x(3)^2)];  
g1 = @(t, x) [0; 0; param.sigma];     

f2 = @(t, x) [x(2); (-x(1) + x(3)); -x(1) - (2*x(2) + 3*x(3))];
g2 = g1;
%start in location 
loc_curr = 1;
x0_curr = x0;
t_curr = 0;

%SDETOOLS for sampling
opts_fw = sdeset('EventsFun', @(t, x) stoch_event_fw(t, x, param),...
    'SDEType', 'Ito',  'DiagonalNoise', 'yes', 'ConstGFUN', 'yes');

opts_bk = sdeset('EventsFun', @(t, x) stoch_event_bk(t, x, param),...
    'SDEType', 'Ito',  'DiagonalNoise', 'yes', 'ConstGFUN', 'yes');

%start the main loop
% out_sim = struct;
% out_sim.sim = {};
% out_sim.jump = {};
osd.locations = {{}, {}};
osd.guards = {{}, {}};

while t_curr <= (Tmax-param.dt) %hack to get the limit right
    time_track_trunc = Tmax - t_curr;
    if loc_curr == 1
        [x_curr,W,TE,YE,WE,IE] = sde_euler(f1,g1,0:param.dt:time_track_trunc,...
        x0_curr,opts_fw);                 
    else
        [x_curr,W,TE,YE,WE,IE] = sde_euler(f2,g2,0:param.dt:time_track_trunc,...
        x0_curr,opts_bk); 
    end
    
    %store data
    time_curr = (1:(size(x_curr, 1)-1))*param.dt;
    time_accum = time_curr + t_curr;
    
    %jump to next location
    t_curr = time_accum(end);
    x0_curr = x_curr(end, :)';
    loc_prev = loc_curr;
    loc_curr = 3 - loc_prev; %toggle between 1 and 2
    
    %document the trajectory and transition
    out_sim = struct;
    out_sim.t = time_accum;
    out_sim.x = x_curr;
    
    jump_curr = struct('t', t_curr, 'x', x0_curr, 'x_jump', x0_curr, 'guard', 3-loc_curr);
    
    osd.locations{loc_prev} = vertcat(osd.locations{loc_prev}, {out_sim});
    osd.guards{loc_prev} = vertcat(osd.guards{loc_prev}, {jump_curr});
    
end

end

