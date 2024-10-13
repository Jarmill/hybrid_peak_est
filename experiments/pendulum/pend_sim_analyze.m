load('pendulum_250.mat');

mm = 0;
for i=1:length(osm)
    curr = osm{i};
    for j = 1:length(curr.sim)
        mm = max(mm, max(max(curr.sim{j}.x(:, 2).^2)));
    end

end