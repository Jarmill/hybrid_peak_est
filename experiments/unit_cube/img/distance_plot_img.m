load('distance_cube_param.mat')
load('traj_deterministic_cube.mat')

%% compute the unsafe set
dist = obj_rec;
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
axl(i).XLim = [-0.6, 1];
end
% view(3)