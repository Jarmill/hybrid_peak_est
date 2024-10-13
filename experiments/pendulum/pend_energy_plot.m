%plot energy

Nth = 90;
Nom = 90;
Om_max = 3;
theta = linspace(0, 2*pi, Nth);
omega = linspace(-Om_max, Om_max, Nom);

[Theta, Omega] = meshgrid(theta, omega);


E_func = @(th, w) 0.5.*w.^2 - cos(th);


E = E_func(Theta, Omega);

%lqr realization
%friction
% b = 0.001;
b = 0;
A = [0 1; 1 -b];
B = [0; -1];

Q = diag([2, 1]);
R = 5;
[K, S, clp] = lqr(A, B, Q, R);

quad_lqr = @(th, w) S(1,1)*(-sin(th)).^2 + 2*S(1,2).*(-sin(th)).*w + S(2, 2)*w.^2;

theta_lqr = linspace(pi/3, 5*pi/3, Nth);
omega_lqr = linspace(-2, 2, Nom);
[Theta_lqr, Omega_lqr] = meshgrid(theta_lqr, omega_lqr);
Quad_lqr = quad_lqr(Theta_lqr, Omega_lqr);
E_lqr = E_func(Theta_lqr, Omega_lqr);


%parameters of hybrid controller
epsilon = 0.1; %decision boundary for swing-up
% delta = 2; %decision boundary for LQR-like control
delta = 1; %decision boundary for LQR-like control

mask_lqr = (E_lqr >= 1-epsilon) &  (E_lqr <= 1+epsilon);
Quad_lqr_mask = Quad_lqr;




% Quad_lqr_mask(~mask_lqr) = NaN;

%plots
% figure(1)
% clf
% surf(Theta, Omega, E);
% xlabel('\theta')
% ylabel('\omega')
% zlabel('E')
% title('pendulum energy landscape')

figure(2)
clf
hold on
contour(Theta, Omega, E, [1-epsilon, 1+epsilon], 'k', 'DisplayName', 'Energy Target')
% contour(Theta, Omega, E, 1+epsilon)
% contourf(Theta, Omega, -E, -[1 1]*[1-epsilon])
% contourf(Theta, Omega, E, [1 1]*[1+epsilon])
% colormap([1.0000    0.6979         0])
% contour(Theta, Omega, -E, -[1 1]*(1-epsilon))
% contour(Theta, Omega, E, [1 1]*(1+epsilon))
contour(Theta_lqr, Omega_lqr, Quad_lqr_mask, linspace(0, delta, 10), 'k', 'DisplayName', 'LQR region')
xlabel('\theta')
ylabel('\omega')
zlabel('E')
title('pendulum energy landscape')