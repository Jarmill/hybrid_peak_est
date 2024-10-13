w_range = [-0.5; 0.5];
c_range = [-1; 1];

    ps_func = @() PM.pend_sample_point_box(c_range, w_range);
    ps = struct('N', 10, 'init', ps_func);
osm = PM.pend_sample_multi(ps, 10);

% x0 = [-1.95489223011601
% 0.269596929272335];
% 
% % x0 = [-1.95489223011601;
% %     0];
% 
% % x0 = [1; 0];
% % x0 = [2; 0.5];
% x0 = [0; 1];
% 
% % x0 = [0; 0.269596929272335];
% 
% osn = PM.pend_sim_nonneg(x0, 10);