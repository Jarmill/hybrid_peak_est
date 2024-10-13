load('det_cube_big.mat')
% load('noisy_cube_big.mat')


    Ru = 0.25;
    Cu = [0.7; 0; 0.5];

dd = @(x) half_sphere_dist(x-Cu, Ru);

figure(2)
clf
hold on
mm = Inf;
for i = 1:length(osd.locations{2})
    curr = osd.locations{2}{i};
    % plot3(curr.x(:, 1), curr.x(:, 2), curr.x(:, 3), 'k');
    % mm = max(mm, max(max(abs(curr.x(:, 1)))));
    % dlist = zeros(length(curr.x), 1)
    for j = 1:length(curr.x)
        mm = min(mm, dd(curr.x(j, :)));
    end
end

mm


%distance 
function dist_out = half_sphere_dist(x_in, R)
    dist_out = half_circ_dist([hypot(x_in(1), x_in(2)), x_in(3)], R);
end

function dist_out = half_circ_dist(x_in, R)
    %return the L2 distance  between the point x_in and the half circle
    %||x_in||^2 <= R^2 intersect x_in(2) <= 0.
    
    if x_in(2) >= 0
        %flat region
        if x_in(1) < -R
            dist_out = hypot(x_in(1)+R, x_in(2))^2;
        elseif x_in(1) > R
            dist_out = hypot(x_in(1)-R, x_in(2))^2;
        else
            dist_out = x_in(2)^2;
        end
    else
        %circle region
        dist_out = max(norm(x_in, 2)-R, 0);
    end

end