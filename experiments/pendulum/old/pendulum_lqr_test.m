%pendulum with unit constants (length, mass, gravity)
rng(300, 'twister')

CONTROLLER = 2;


%lqr realization
A = [0 1; 1 0];
B = [0; -1];

% Q = diag([2, 1]);
Q = diag([2, 1]);
R = 5;
%  R = 1;

[K, S, clp] = lqr(A, B, Q, R);


%this is a stable choice when |th-pi|<=pi/4 and |w| <= 1.
%u_lqr = @(x) K*[sin(x(1) - pi); x(2)];
u_lqr = @(x) K*[-sin(x(1)); x(2)];


%swing-up
% k = 1; %swing-up gain
% k = 0.5; %swing-up gain
k = 0.75; %swing-up gain
u_swing = @(x) -k * x(2) * ((0.5*x(2)^2 - cos(x(1))) - 1);



if CONTROLLER == 0
    %lqr    
    th_max = pi/4;
    w_max = 1;
    u = u_lqr;
elseif CONTROLLER == 1    
    th_max = pi;
    w_max = 1;

    u = u_swing;
else
    th_max = pi;
    w_max = 1.5;

    eps = 0.1;        
    delta = 1;
%     delta = 2;
    u = @(x) u_hybrid(x, S, delta, eps, u_lqr, u_swing);
end

%% Simulate
% Nsample = 300;
Tmax_sim = 20;
Nsample = 100;
sampler = @()  [1e-2; 0]; 
% sampler = @() [pi + (2*rand() - 1)* th_max; (2*rand()-1)*w_max];
%out_sim = pend_sampler(sampler, Nsample, Tmax_sim);
out_sim = pend_sampler(sampler, Nsample, Tmax_sim, [], 0, u);

%% Plot
figure(1)
clf
for i = 1:length(out_sim)
    subplot(5, 1, 1)
    if i== 1
        hold on
        title('Cosine of Angle')
        xlabel('time')
        ylabel('cos(\theta)')
    end
    plot(out_sim{i}.t, out_sim{i}.x(:, 1), 'c')

    subplot(5, 1,2)
    if i== 1
        hold on
        title('Sine of Angle')
        xlabel('time')
        ylabel('sin(\theta)')
    end
    plot(out_sim{i}.t, out_sim{i}.x(:, 2), 'c')

    subplot(5, 1,3)
    if i== 1
        hold on
        title('Angular Velocity')
        xlabel('time')
        ylabel('\omega')
    end
    hold on
    plot(out_sim{i}.t, out_sim{i}.x(:, 3), 'c')
    
    out_sim{i}.energy = 0.5*out_sim{i}.x(:, 3).^2 - out_sim{i}.x(:, 1);
    
    subplot(5, 1,4)
    if i== 1
        hold on
        title('Energy')
        xlabel('time')
        ylabel('E')
    end
    hold on
    plot(out_sim{i}.t, out_sim{i}.energy, 'c')
    
    
    out_sim{i}.u = zeros(size(out_sim{i}.t));
    for k = 1:length(out_sim{i}.u)
        out_sim{i}.u(k) = u(trig_process(out_sim{i}.x(k, :)));
    end
    subplot(5, 1,5)
    if i== 1
        hold on
        title('Input')
        xlabel('time')
        ylabel('u')
    end
    hold on
    plot(out_sim{i}.t, out_sim{i}.u, 'c')
end

subplot(5, 1, 1)
plot(xlim, [-1, -1], ':k')
subplot(5, 1, 2)
plot(xlim, [0, 0], ':k')
subplot(5, 1, 3)
plot(xlim, [0, 0], ':k')
subplot(5, 1, 4)
plot(xlim, [1, 1], ':k')
subplot(5, 1, 5)
plot(xlim, [0, 0], ':k')

figure(2)
clf
hold on
for i = 1:length(out_sim)
%     plot(out_sim{i}.theta, out_sim{i}.x(:, 3), 'c')
%     plot3(out_sim{i}.t, out_sim{i}.theta, out_sim{i}.x(:, 3), 'c')
%     scatter3(0, out_sim{i}.theta(1), out_sim{i}.x(1, 3), 200, 'ok')
plot3(out_sim{i}.x(:, 1), out_sim{i}.x(:, 2), out_sim{i}.x(:, 3), 'c')
scatter3(out_sim{i}.x(1, 1), out_sim{i}.x(1, 2), out_sim{i}.x(1, 3), 50, 'ok')
end
scatter3(-1, 0, 0, 400, 'rs', 'LineWidth', 2)
hold off
xlabel('cos(\theta)')
ylabel('sin(\theta)')
zlabel('omega(\theta)')
% xlim([0, 2*pi])
view(3)
title('Swing-up Pendulum Trajectories', 'fontsize', 16)


figure(3)
clf
hold on
for i = 1:length(out_sim)
    theta_nan = out_sim{i}.theta;
    theta_nan_diff = [0; abs(diff(theta_nan))];
    theta_nan(theta_nan_diff > 5.5) = NaN;
%     plot(out_sim{i}.theta, out_sim{i}.x(:, 3), '-c')
    plot(theta_nan, out_sim{i}.x(:, 3), '-c')
%     plot3(out_sim{i}.t, out_sim{i}.theta, out_sim{i}.x(:, 3), 'c')
    scatter(out_sim{i}.theta(1), out_sim{i}.x(1, 3), 20, 'ok')
end

%energy contour
E_gap_func = @(th, w) 0.5.*w.^2 - cos(th) - 1;
quad_func = @(th, w) S(1,1).*(-sin(th)).^2 + 2.*S(1,2).*(-sin(th)).*w + S(2, 2).*w.^2;

% quad_func_mask = @(th, w) quad_func(th, w) .* (E_gap_func(th, w) <= eps) .* (E_gap_func(th, w) >= -eps);
fcontour(E_gap_func, [0, 2*pi, -2.5, 2.5], 'k', 'levellist', [-eps, eps])


%LQR control contour
% theta_lqr = linspace(0, 2*pi, 100);
% omega_lqr = linspace(-2.5, 2.5, 100);
% [Theta_lqr, Omega_lqr] = meshgrid(theta_lqr, omega_lqr);
% Quad_lqr = quad_func(Theta_lqr, Omega_lqr);
% E_lqr = E_gap_func(Theta_lqr, Omega_lqr);
% mask_lqr = (E_lqr >= -eps) &  (E_lqr <= eps);
% Quad_lqr_mask = Quad_lqr;
% contour(Theta_lqr, Omega_lqr, Quad_lqr_mask, delta, '--k', 'DisplayName', 'LQR region')


fcontour(quad_func, [pi/2, 3*pi/2, -2.5, 2.5], '--k', 'levellist', [delta])
% fcontour(quad_func, [0, 2*pi, -2.5, 2.5], '--k', 'levellist', delta)
% fcontour(quad_func_mask, [0, 2*pi, -2.5, 2.5], 'k', 'levellist', delta)
scatter(pi, 0, 300, 'rs', 'linewidth', 1.5)
xlabel('\theta')
ylabel('\omega')
title('Swing-up Pendulum Trajectories', 'fontsize', 16)


function x = trig_process(x_trig)
    x = [atan2(x_trig(2), x_trig(1)); x_trig(3)];
end
    
function u_out = u_hybrid(x, S, delta, eps, u_lqr, u_swing)
%hybrid controller.
%if abs(energy -1) >= eps, use u_swing
%if abs(energy-1) <= eps and norm(x)^2 <= delta, use u_lqr
%else, wait
energy_gap = 0.5.*x(2).^2 - cos(x(1)) - 1;

if abs(energy_gap) >= eps
    u_out = u_swing(x);
else
    %if norm(x)^2 <= delta
    %sx = [sin(x(1) - pi); x(2)];
    sx = [-sin(x(1)); x(2)];
    Sform = sx'*S*sx;
    if Sform <= delta
        u_out = u_lqr(x);
    else
        u_out = 0;
    end
end


%three states:
%swing-up
%wait
%lqr

end