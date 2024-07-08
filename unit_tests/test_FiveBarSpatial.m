rob = FiveBarSpatial();

% diagrams.setup; hold on;
% rob.plot(zeros(6,1), zeros(6,1));
% diagrams.redraw; hold off;

%%

% q = rand_angle([6 1])
% [R_0T, p_0T] = fwdkin(rob.kinA, q);

R_0T = eye(3);
p_0T = 0.3*[1;1;1];

[Q_A, Q_B, is_LS_A, is_LS_B] = rob.IK(R_0T, p_0T)

%%
i = 5;
diagrams.setup; hold on;
rob.plot(Q_A(:,i), Q_B(:,i));
diagrams.arrow([0;0;0], p_0T);
diagrams.redraw; hold off;