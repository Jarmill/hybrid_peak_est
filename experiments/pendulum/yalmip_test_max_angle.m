S = [   12.6422   10.9161
   10.9161   10.6846];

x = sdpvar(2, 1);

cons = [x'*S*x <= 1];
sol = optimize(cons, x(1), sdpsettings('solver', 'mosek', 'verbose', 0));
x_rec = double(x);
sin_theta = -x_rec(1);
theta_rec = asin(sin_theta)