%find the maximum cosine of angle in LQR control range
xs = sdpvar(3, 1); %[sin(theta), w]
% xs_lim = [1; ang_max];
L = chol(PM.S);
F = [norm(xs(1:2), 2) <= 1;
     norm(L*[-xs(2); xs(3)], 2) <= sqrt(PM.delta);
     xs(3) <=  ang_max; 
     xs(3) >= -ang_max];
h = -xs(2);

sol = optimize(F, h, sdpsettings('solver', 'mosek', 'verbose', 0));

double(xs)

% sx = double(xs(1));
% cx = -sqrt(1-sx^2)
% theta_lqr = atan2(sx, cx)*180/pi