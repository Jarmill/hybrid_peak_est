%pendulum with unit constants (length, mass, gravity)

%find maximum angular velocity with hybrid controller
mset clear
rng(300, 'twister')

%states/variables
mpol('t')
mpol('x', 3, 1);
c = x(1);
s = x(2);
w = x(3);


%objective to maximize
objective = x(3)^2;        %maximum angular velocity
%objective = -c;
%objective = 0.5*w^2 - c;    %maximum energy


%support
w_max = 2.5;
T = 10;
Xsupp = [c^2 + s^2 == 1; w^2 <= w_max^2];
% Xsupp = [c^2 + s^2 <= 1; c^2 + s^2 >= 1];
%% Formulate dynamics
%friction
b= 0;

energy = 0.5*w^2 - c;
energy_gap = energy - 1;
f_wait = [-s*w; c*w; -s - b*w];

%lqr realization
A = [0 1; 1 -b];
B = [0; -1];

Q = diag([2, 1]);
R = 5;

[K, S, clp] = lqr(A, B, Q, R);

%this is a stable choice when |th-pi|<=pi/4 and |w| <= 1.
u_lqr = @(x) K*[-sin(x(1)); x(2)];

u_lqr_poly = K*[-s; w];
%lqr decision
quad_lqr = [-s; w]'*S*[-s; w];

f_lqr = f_wait + [0; 0; 1]*u_lqr_poly;

%swing-up
k = 1; %swing-up gain
% k = 1.2; %swing-up gain

u_swing = @(x) -k*x(2)*(0.5*x(2)^2 - cos(x(1)));
u_swing_poly = -k * w * energy_gap;
f_swing = f_wait + [0; 0; 1]*u_swing_poly;

%parameters of hybrid controller
epsilon = 0.1; %decision boundary for swing-up
delta = 1; %decision boundary for LQR-like control

% X_base = [c^2 + s^2 == 1; w^2 <= w_max^2]
X_up = (energy <= (1 - epsilon));
X_down = (energy >= (1+epsilon));
% X_swing_pos = energy_gap >= epsilon;
% X_swing_neg = energy_gap <= -epsilon;
X_wait  = [energy <= (1+epsilon); energy >= (1-epsilon); quad_lqr >= delta];
X_lqr   = [energy <= (1+epsilon); energy >= (1-epsilon); quad_lqr <= delta];
% X_wait  = [energy_gap <= epsilon; energy_gap >= -epsilon; quad_lqr >= delta];
% X_lqr   = [energy_gap <= epsilon; energy_gap >= -epsilon; quad_lqr <= delta];


%wrap up the dynamics

f = {T*f_swing, T*f_swing, T*f_wait, T*f_lqr};
X = {X_up, X_down, X_wait, X_lqr};

% th_max = pi-0.01;
% w_max = 1.5;

w_max_0 = 0.75;
X0 = [s^2 + c^2 == 1; w^2 <= w_max_0^2];

%formulate problem
p_opt = peak_options;
p_opt.var.t = t;
p_opt.var.x = x;

p_opt.dynamics = struct;
p_opt.dynamics.f = f;
p_opt.dynamics.X = X;

p_opt.state_init = X0;
p_opt.state_supp = Xsupp;
p_opt.Tmax = 1;
% p_opt.box = [-1, 1; -1, 1; -4, 4];
p_opt.box = [];
p_opt.scale = 0;

p_opt.obj = objective;
%p_opt.Tmax = 20;
% order = 2;
% out = peak_estimate(p_opt, order);


 
%% Simulate
%Nsample = 20;
Nsample = 100;
Tmax_sim = 20;

u = @(x) u_hybrid(x, S, epsilon, delta, u_lqr, u_swing);

HALF_ANGLE = 0;
% if HALF_ANGLE
%     sampler = @() [pi - rand()* th_max; (2*rand()-1)*w_max];
% else
%     sampler = @() [pi + (2*rand() - 1)* th_max; (2*rand()-1)*w_max];
% end
sampler = @()[2*pi*rand(); (2*rand()-1)*w_max];
out_sim = pend_sampler(sampler, Nsample, Tmax_sim, out.dynamics.nonneg, b, u);

if (out.optimal == 1)
    x0_angle = [atan2(out.x0(2), out.x0(1)); out.x0(3)];
    out_sim_peak = pend_sampler(x0_angle, 1, Tmax_sim, out.dynamics.nonneg, b, u);
    nplot = nonneg_plot(out_sim, out_sim_peak);
    pplot = pend_plot_3(out_sim, out_sim_peak);
%     splot = state_plot_3(out, out_sim, out_sim_peak);
% %     splot = state_plot_N(out, out_sim, out_sim_peak);
else
    nplot = nonneg_plot(out, out_sim);
%     splot = state_plot_2(out, out_sim);
%     splot3 = state_plot_3(out, out_sim);
    %splot = state_plot_N(out, out_sim);
    pplot = pend_plot_3(out, out_sim);
end

function u_out = u_hybrid(x, S, eps, delta, u_lqr, u_swing)
%hybrid controller.
%if abs(energy -1) >= eps, use u_swing
%if abs(energy-1) <= eps and norm(x)^2 <= delta, use u_lqr
%else, wait
energy_gap = 0.5*x(2)^2 - cos(x(1)) - 1;

if abs(energy_gap) >= eps
    u_out = u_swing(x);
else
    %if norm(x)^2 <= delta
%     sx = [sin(x(1) - pi); x(2)];
    sx = [-sin(x(1)); x(2)];
    Sform = sx'*S*sx;
    if Sform <= delta
        u_out = u_lqr(x);
    else
        u_out = 0;
    end
end

end