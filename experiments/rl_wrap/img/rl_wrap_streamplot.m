R = 1;
N = 31;
[x, y] = meshgrid(linspace(-R, R, N));

%original
% k = 0.5;
k = 0.6;
dx = (-y + x.*y+k);
dy = (-y  - x + x.^3);




figure(444)
clf

% quiver(x, y, dx, dy)
% verts = stream2(x, y, dx, dy, x, y);
% streamline(verts)
streamslice(x, y, dx, dy)
axis square