%% variables and dynamics

x = sdpvar(3, 1);       %state
t = sdpvar(1, 1);       %time
w = [];
gamma = sdpvar(1, 1);   %bound

alpha_fw = sdpvar(1, 1);
alpha_bk = sdpvar(1, 1);

% Tmax = 5;
Tmax = 20;
order = 2;
d = 2*order;


%dynamics
f1 = [x(2); -x(1) + x(3); x(1) + (2*x(2) + 3*x(3))*(1 + x(3)^2)];
f2 = [x(2); -x(1) + x(3); -x(1) - 2*x(2) - 3*x(3) ];

%peak maximizing objectives
% p1 = x(1)^2;
p1 = [];
p2 = x(1)^2;


%% support sets

r2 = sum(x.^2);
% K = 1;
critfw = x(1)^2 + x(2)^2 + x(3)^2;

% L = [6; 12; 12];
% L = [5.5; 4; 10.5];
box_X = 1.5.^2-x.^2;

    R0 = 0.2;
    R1 = 1;

X01 = struct('ineq', [], 'eq', [r2-R0^2]);
X02 = struct('ineq', [], 'eq', []);
% X1 = struct('ineq', [box_X; 1.01-critfw; t*(1-t)], 'eq', []);
X1 = struct('ineq', [R1^2-critfw; t*(1-t)], 'eq', []);
X2 = struct('ineq', [box_X; r2 - R0^2; t*(1-t)], 'eq', []);
X1_lie = struct('ineq', X1.ineq, 'eq', []);
X2_lie = struct('ineq', X2.ineq, 'eq', []);

%guards
% Xfw = struct('ineq', [box_X; 1.01-critfw; critfw - 0.99; t*(1-t)], 'eq', []);
Xfw = struct('ineq', [t*(1-t)], 'eq', [critfw-R1^2]);

critbk = r2;
Xbk = struct('ineq', [t*(1-t)], 'eq', [r2-R0^2]);

%% polynomials
%symmetry (x, w) <-> (-x, -w)
% 
% var_count_init = [1, 0, 2, 0]; %x
% var_count_std = [0, 0, 3, 1]; %x t
% var_count_lie = [0, 0, 3, 1]; %x t
% var_count_guard = [1, 0, 2, 1]; %x t

% mom_out = sym_moments(var_count_std, order);
% v_monom = recovermonoms(mom_out.monom_int,[x;t]);
% cv1 = sdpvar(length(v_monom), 1);
% cv2 = sdpvar(length(v_monom), 1);
% 
% v1 = cv1'*v_monom;
% v2 = cv2'*v_monom;

v1 = polynomial([t;x], d);
v2 = polynomial([t;x], d);

%% psatz constraints
% epsilon = 1e-8;
epsilon = 0;
%initial
v1_init = replace(v1, t, 0);
[p_init, cons_init, coeff_init] = constraint_psatz(gamma - v1_init-epsilon, X01, x, d);

%peak
% [p_cost1, cons_cost1, coeff_cost1] = constraint_psatz(v1 - p1-epsilon, X1, [x; t], order);
[p_cost2, cons_cost2, coeff_cost2] = constraint_psatz(v2 - p2-epsilon, X2, [x; t], d);

%lie derivative
Lv1 = jacobian(v1, x)*f1*Tmax+ jacobian(v1, t);
Lv2 = jacobian(v2, x)*f2*Tmax+ jacobian(v2, t);
[p_lie1, cons_lie1, coeff_lie1] = constraint_psatz(-Lv1-epsilon, X1_lie, [x; t], d+2);
[p_lie2, cons_lie2, coeff_lie2] = constraint_psatz(-Lv2-epsilon, X2_lie, [x; t], d+2);

%guards
[p_fw, cons_fw, coeff_fw] = constraint_psatz(v1 - v2 + alpha_fw-epsilon, Xfw, [x; t], d);
[p_bk, cons_bk, coeff_bk] = constraint_psatz(v2 - v1 + alpha_bk-epsilon, Xbk, [x; t], d);

%zeno? 
ZENO_CAP = true;
N_fw = 5;
N_bk = 5;

if ZENO_CAP
    cons_zeno = [alpha_fw >= 0; alpha_bk >= 0];
else
    cons_zeno = [alpha_fw == 0; alpha_bk == 0];
end
    

%accumulate the constraints
cons = [cons_init:'Init'; cons_cost2:'Cost 2'; ... %cons_cost1:'Cost 1'; 
    cons_lie1:'Lie 1'; cons_lie2:'Lie 2'; cons_fw:'Guard 1->2'; cons_bk:'Guard 2->1';
    cons_zeno:'Zeno'];

coeff = [coeff_fw; coeff_bk; coeff_lie1; coeff_lie2; coeff_cost2; coeff_init;...    
    gamma; alpha_fw; alpha_bk];

%% solve the system



objective = gamma + N_fw * alpha_fw + N_bk * alpha_bk;


MP = false;

if MP
%     opts = sdpsettings();
opts = sdpsettings('solver','sdpa_gmp',...
    'sdpa_gmp.epsilonStar', 10^(-25), ...
    'sdpa_gmp.epsilonDash', 10^(-25), ...
    'sdpa_gmp.lambdaStar', 10^(4), ...
    'sdpa_gmp.betaStar', 0.1, ...
    'sdpa_gmp.betaBar',0.3, ...
    'sdpa_gmp.gammaStar',0.7,  ...
    'sdpa_gmp.lowerBound',0,  ... %changed lower bound to 0
    'sdpa_gmp.upperBound',10^25,  ...
    'sdpa_gmp.omegaStar',2,  ...
    'sdpa_gmp.maxIteration',200, ... 
    'sdpa_gmp.precision',250);
else
     opts = sdpsettings('solver', 'mosek','savesolverinput', 1);
end

[sol, monom, Gram, residual] = solvesos(cons, objective, opts, [coeff]);

% sol = optimize(cons, objective, opts);

% value(objective)
fprintf('abs(x1) bound: %0.4f \n', sqrt(value(objective)))


%% save
v1_rec = value(cv1)'*v_monom;
v2_rec = value(cv2)'*v_monom;
v1_func = polyval_func(v1_rec, [x; t]);
v2_func = polyval_func(v2_rec, [x; t]);
alpha_fw_rec = value(alpha_fw);
alpha_bk_rec = value(alpha_bk);
gamma_rec = value(gamma);

alpha_rec = [alpha_fw_rec; alpha_bk_rec];

Lv1_rec = jacobian(v1_rec, x)*f1*Tmax+ jacobian(v2_rec, t);
Lv2_rec = jacobian(v2_rec, x)*f2*Tmax + jacobian(v2_rec, t);

f1_func = polyval_func(f1, [w; x; t]);
f2_func = polyval_func(f2, [w; x; t]);

Lv1_func = polyval_func(Lv1, [w; x; t]);
Lv2_func = polyval_func(Lv2, [w; x; t]);

nn_1_rec = [gamma_rec - v1_rec;
          v1_rec - p1;
           -Lv1_rec];
          
nn_2_rec = [v2_rec - p2; 
            -Lv2_rec];
      
nn_fw = [v1_rec - v2_rec + alpha_fw_rec];
nn_bk = [v2_rec - v1_rec + alpha_bk_rec];

nn_1_func = polyval_func(nn_1_rec, [w; x; t]);
nn_2_func = polyval_func(nn_2_rec, [w; x; t]);
nn_fw_func = polyval_func(nn_fw, [x; t]);
nn_bk_func = polyval_func(nn_bk, [x; t]);

obj_rec = value(objective);
% save('dual_deterministic_5.mat', 'obj_rec', 'nn_fw_func', 'nn_bk_func', ...
% 'nn_1_func', 'nn_2_func', 'v1_func', 'v2_func', 'Lv1_func', 'Lv2_func',...
% 'f1_func', 'f2_func', 'Lv1_func', 'Lv2_func', 'gamma_rec', 'alpha_rec', 'sol', ...
% 'order', 'Tmax');

% f1f = polyval_func(f1, [w;x]);