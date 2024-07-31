% Test subproblem_8th_order.dual_quadratic

% Generate random coefficients (except constants) and random x_1, x_2

P_a = rand(1, 3);
P_b = rand(1, 3);
P_c = rand(1, 3);
P_c(3) = NaN;
Q_a = rand(1, 3);
Q_b = rand(1, 3);
Q_c = rand(1, 3);
Q_c(3) = NaN;
x_1 = rand;
x_2 = rand;

P_a_i = P_a(1)*x_1^2 + P_a(2)*x_1 + P_a(3);
P_b_i = P_b(1)*x_1^2 + P_b(2)*x_1 + P_b(3);
P_c_i_tilde = P_c(1)*x_1^2 + P_c(2)*x_1;
Q_a_i = Q_a(1)*x_1^2 + Q_a(2)*x_1 + Q_a(3);
Q_b_i = Q_b(1)*x_1^2 + Q_b(2)*x_1 + Q_b(3);
Q_c_i_tilde = Q_c(1)*x_1^2 + Q_c(2)*x_1;

% Find the constants Pc_3 and Qc_3
P_c(3) = -(P_a_i*x_2^2 + P_b_i*x_2 + P_c_i_tilde);
Q_c(3) = -(Q_a_i*x_2^2 + Q_b_i*x_2 + Q_c_i_tilde);

% Double check two equations equal 0
% Pa(x1) x2^2 + Pb(x1) x2 + Pc(x1) = 0
% Qa(x1) x2^2 + Qb(x1) x2 + Qc(x1) = 0
P_c_i = P_c(1)*x_1^2 + P_c(2)*x_1 + P_c(3);
Q_c_i = Q_c(1)*x_1^2 + Q_c(2)*x_1 + Q_c(3);

assert(abs(P_a_i*x_2^2 + P_b_i*x_2 + P_c_i) < 1e-9);
assert(abs(Q_a_i*x_2^2 + Q_b_i*x_2 + Q_c_i) < 1e-9);

% Solve the system of equations
[x1_t, x2_t] = subproblem_8th_order.dual_quadratic(P_a, P_b, P_c, Q_a, Q_b, Q_c);
[x_1 x_2]
[x1_t x2_t]