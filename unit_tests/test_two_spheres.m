% Test two_spheres.m

% Generate random p1 p2 p3 p4 R1 R2
p1 = rand_vec;
p2 = rand_vec;
p3 = rand_vec;
p4 = rand_vec;
k1 = rand_normal_vec;
k2 = rand_normal_vec;
theta1 = rand_angle;
theta2 = rand_angle;
R1 = rot(k1, theta1);
R2 = rot(k2, theta2);

% Calculate d1, d2
d1 = norm(R1*p1 - R2*p2);
d2 = norm(R1*p3 - R2*p4);

[t1, t2] = subproblem_8th_order.two_spheres(p1, p2, p3, p4, k1, k2, d1, d2);
[theta1 theta2]
[t1 t2]

% test 1 soln
norm(rot(k1, t1(1)) * p1 - rot(k2, t2(1)) * p2) - d1
norm(rot(k1, t1(1)) * p3 - rot(k2, t2(1)) * p4) - d2