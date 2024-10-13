mset clear

mpol('t')
mpol('x', 3, 1)

vars.t = t;
vars.x = x;

supp = [t*(1-t) >= 0; x.^2 <= 1];
supp0 = [t == 0; x == 0.5];

% mu = meas_base(1, vars, supp);

loc = location(1, vars, supp, supp0);

g = guard(1, vars, loc, loc, supp, 3*x);

[ms, md] = g.liou_reset(3)
% f = -vars.x + 0.4*vars.x(1);
% 
% Liou = loc.liou_con(4, f);