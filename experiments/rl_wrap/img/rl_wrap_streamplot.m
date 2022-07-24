R = 1;
N = 31;
[x, y] = meshgrid(linspace(-R, R, N));

%original
% dx = (-y + x.*y);
% dy = (-y  - x + x.^3);

dx = (-y + x.*y+0.5);
dy = (-y  - x + x.^3);


figure(444)
clf

% quiver(x, y, dx, dy)
% verts = stream2(x, y, dx, dy, x, y);
% streamline(verts)
streamslice(x, y, dx, dy)
axis square