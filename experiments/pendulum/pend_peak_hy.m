%find maximum angular velocity with hybrid controller
mset clear
rng(300, 'twister')

%states/variables
mpol('t', 1, 1);
mpol('x', 3, 1);
c = x(1);
s = x(2);
w = x(3);
vars = struct('t', t, 'x', x);

%objective to maximize
objective = x(3)^2;        %maximum angular velocity

E = 0.5*w^2 - c;
f_wait = [-s*w; c*w; -s];

%lqr realization
A = [0 1; 1 0];
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

%swing up control
% k = 1;
k = 1.5;
u_swing_poly = -k * w * (E-1);
f_swing = f_wait + [0; 0; 1]*u_swing_poly;

%parameters of hybrid controller
epsilon = 0.1; %decision boundary for swing-up
delta = 1; %decision boundary for LQR-like control

Tmax = 10; %maximum time
% Tmax = 5;

%% Location Definitions

%initial point 
%[t; c; s; w]
x0 = [0; cos(1e-3); sin(1e-3); 0];

w_max = 2.5; %prior assumption of maximum angular velocity

X = [t*(1-t) <= 1; c^2+s^2 == 1; w^2 <= w_max^2];
% X = [t*(1-t) <= 1; c^2+s^2 <= 1; c^2 + s^2 >= 1; w^2 <= w_max^2];
X_pump_up = [X; E <= 1-epsilon];   %add energy
X_pump_down = [X; E >= 1+epsilon]; %reduce energy
X_wait = [X; E <= 1+epsilon; E >= 1- epsilon; quad_lqr >= delta];
X_lqr = [X; E <= 1+epsilon; E >= 1- epsilon; quad_lqr <= delta];

loc_up   = location(vars, X_pump_up,   x0, f_swing, objective, 1);
loc_wait = location(vars, X_wait, [], f_wait, objective, 2);
loc_lqr  = location(vars, X_lqr,  [], f_lqr, objective, 3);
loc_down = location(vars, X_pump_down, [], f_swing, objective, 4);
locations = {loc_up; loc_wait; loc_lqr; loc_down};
% locations = {loc_up; loc_wait; loc_lqr};
%% Guard Definition

%identity reset map
Reset = x;

X_up_wait = [X; E == 1 + epsilon; quad_lqr >= delta];
X_up_lqr = [X; E == 1 + epsilon; quad_lqr <= delta];
X_down_wait = [X; E == 1 - epsilon; quad_lqr >= delta];
X_down_lqr = [X; E == 1 - epsilon; quad_lqr <= delta];
X_wait_lqr = [X; E <= 1 + epsilon; E >= 1 - epsilon; quad_lqr == delta];


g_up_wait = guard(1, vars, loc_up, loc_wait, X_up_wait, Reset);
g_up_lqr = guard(2, vars, loc_up, loc_lqr, X_up_lqr, Reset);
g_down_wait = guard(3, vars, loc_down, loc_wait, X_down_wait, Reset);
g_down_lqr = guard(4, vars, loc_down, loc_lqr, X_down_lqr, Reset);
g_wait_lqr = guard(5, vars, loc_wait, loc_lqr, X_wait_lqr, Reset);

guards = {g_up_wait; g_up_lqr; g_down_wait; g_down_lqr; g_wait_lqr};
% guards = {g_up_wait;  g_up_lqr; g_wait_lqr};
%% Peak Manager
PM = peak_manager_hy(locations, guards);

%run peak
order = 3;
% [objective_mom, mom_con, supp_con] =  PM.peak_cons(order);
sol = PM.peak(order);
w_max_est = sqrt(sol.obj_rec)