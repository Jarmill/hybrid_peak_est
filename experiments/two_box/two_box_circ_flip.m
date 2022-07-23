rng(25, 'twister')
mset clear
mpol('t');
mpol('x', 2, 1);
vars = struct('t', t, 'x', x);

% X = [t*(1-t) <= 1; x.^2 <= 1];

Tmax = 3;

%box 1
% X01 = [0; 0.3; 0.5];
X01 = [t == 0;  (x - [0.3; 0.3]).^2 <= 0.01];
% f1 =  x;
f1 = [x(1) + 3*x(2); -3*x(1) + x(2) + x(2)^3 - x(1)*x(2)^2];
% p1 = x(2);
p1 = [];
lsupp1 = loc_support(vars);
lsupp1 = lsupp1.set_box(1);
lsupp1.X_init = X01;
lsupp1.Tmax = Tmax;
loc1 = location(lsupp1, f1, p1, 1);

%box 2
% X02 = [0; 0.5; 0.5];
% X02 = [];
X02 = [t == 0;  (x - [-0.4; 0]).^2 <= 0.01];
f2 = -[x(1) + x(2) + 0.3*x(2)^3; -3*x(1) + x(2)]; % - x(2)^3
% p2 = -x;
p2 = x(1);
% p2 = [];
lsupp2 = loc_support(vars);
lsupp2 = lsupp2.set_box(1);
lsupp2.Tmax = Tmax;
loc2 = location(lsupp2, f2, p2, 2);

%guard

guard_mom_sub = 0;

%need to figure out how to index the dual variables when there are equality
%support constraints 
% Xgvert = [t*(1-t) <= 1; x(1)^2 - 1 == 0; x(2)^2 <= 1];
% Xghor  = [t*(1-t) <= 1; x(2)^2 - 1 == 0 ;  x(1)^2 <= 1];
Xgtop = [t*(1-t) <= 1; x(2)==1 ;  x(1)^2 <= 1];
Xgbot = [t*(1-t) <= 1; x(2)==-1 ;  x(1)^2 <= 1];
Xgright = [t*(1-t) <= 1; x(1)==1 ;  x(2)^2 <= 1];
Xgleft = [t*(1-t) <= 1; x(1)==-1 ;  x(2)^2 <= 1];
 

Rgcirc = 0.4;
Xgcirc = [t*(1-t) <= 1; Rgcirc^2 - sum(x.^2) >= 0];

gtop =    guard(1, vars, loc1, loc2, Xgtop, x);
gright =  guard(2, vars, loc1, loc2, Xgright, [1;-1].*x);
gbot =  guard(3, vars, loc1, loc2, Xgbot, x);
gleft = guard(4, vars, loc1, loc2, Xgleft, [1;-1].*x);

% gcirc = guard(5, vars, loc2, loc1, Xgcirc, [1;-1].*x);


% gcirc


%run peak
PM =  peak_manager_hy({loc1, loc2}, {gtop, gright, gbot, gleft});

% order = 7;
order = 4;
[sol, PM] = PM.run(order);
sol
% os = PM.sample_traj(0, [0.3; 0.5], 1, 1);

