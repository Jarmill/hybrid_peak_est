% load('distance_cube_param.mat')
load('distance_cube_small_L_param.mat')
load('traj_deterministic_cube.mat')

%% compute the unsafe set
% dist = obj_rec;
% dist = sqrt(0.051026222539867);

%near parameters
% dist = sqrt(0.001040262107855);
% dist = 0.05;

%close but possibly numerically unstable
dist = sqrt(0.005118592284530);
% Cu = [-0.5; -0.3; 0.5];
Cu = [-0.5; -0.5; 0.5];
Ru = 0.4;

%further away, but with a safety-margin certificate of safety
dist = sqrt(0.007941915368103);
Ru = 0.25;
Cu =  [-0.6; -0.5; 0.5];

%quarter-torus
R = Ru;  %radius of half-sphere
r = dist; %inner radius of torus

Nth = 30;
Nphi = 40;

% phi = linspace(4*pi/2, 3*pi, Nphi);
phi = linspace(0, 2*pi, Nphi);  %perimiter
% theta = linspace(0, 2*pi, Nth); %azimuth
theta = linspace(0, -pi/2, Nth);

[Phi, Theta] = meshgrid(phi, theta);

%
X_torus = (R+r.*cos(Theta)).*cos(Phi);
Y_torus = (R+r.*cos(Theta)).*sin(Phi);
Z_torus = r.*sin(Theta);

X_sphere = (R+r)*cos(Theta).*cos(Phi);
Y_sphere = (R+r)*cos(Theta).*sin(Phi);
Z_sphere = -(R+r)*sin(Theta);

rscale = R/(R+r);
circ = R*[cos(phi); sin(phi)];
z_circ = -r * ones(size(phi));
% [X_torus, Y_torus, Z_torus] = pol2cart(Theta, Phi, Z_top);


%% set up the cube plotter
CPlot = cube_plotter(osm, osd, R0, R1);
[F, ax1, ax2] = CPlot.cube_plot();

%% plot the sets
% figure(2)
% clf

axl = [ax1, ax2];
falpha = 0.4;
for i = 1:length(axl)
hold on



%unsafe set
surf(axl(i), rscale*X_sphere + Cu(1), rscale*Y_sphere + Cu(2), rscale*Z_sphere+ Cu(3), 'FaceColor', 'r');
patch(axl(i), 'XData', circ(1,:)+ Cu(1),  'YData', circ(2, :)+ Cu(2),  'ZData', zeros(size(phi))+ Cu(3), 'FaceColor', 'r');

%distance contour

surf(axl(i), X_torus+ Cu(1), Y_torus+ Cu(2), Z_torus+ Cu(3), 'FaceColor', 'r', 'FaceAlpha', falpha, 'edgecolor', 'none');
surf(axl(i), X_sphere+ Cu(1), Y_sphere+ Cu(2), Z_sphere+ Cu(3), 'FaceColor', 'r', 'FaceAlpha', falpha, 'edgecolor', 'none');

patch(axl(i), 'XData', circ(1,:)+ Cu(1), 'YData', circ(2, :)+ Cu(2), 'ZData', z_circ+ Cu(3), 'FaceColor', 'r', 'FaceAlpha', falpha, 'edgecolor', 'none');
axl(i).XLim = [-1, 0.6];
axl(i).YLim = [-1, 0.6];



end

% %optimal point
% scatter3(ax1, mm{1}.x0(1), mm{1}.x0(2),mm{1}.x0(3), 300, 'blue', 'o', 'LineWidth', 2);
% scatter3(ax1, mm{1}.xp(1), mm{1}.xp(2),mm{1}.xp(3), 300, 'blue', '*', 'LineWidth', 2);
% scatter3(ax1, mm{1}.y(1), mm{1}.y(2),mm{1}.y(3), 300, 'blue', 's', 'LineWidth', 2);
% % view(3)