rob = Hexapod();

% q0 = zeros(6,1);

% diagrams.setup; hold on;
% rob.plot(q0, q0, q0, q0, q0, q0);
% diagrams.redraw; hold off;

[Q_A, Q_B, Q_C, Q_D, Q_E, Q_F] = rob.IK(eye(3), [0;0;2]);
i = 1;

diagrams.setup; hold on;
rob.plot(Q_A(:,i), Q_B(:,i), Q_C(:,i), Q_D(:,i), Q_E(:,i), Q_F(:,i));
diagrams.arrow([0;0;0], [0;0;2]);
diagrams.redraw; hold off;

%%
[R_t, T_t] = fwdkin(rob.kinC, Q_C(:,i))

%%

q0 = zeros(6,1);


diagrams.setup; hold on;
rob.plot([pi/2 0 2 0 0 0]', q0, q0, q0, q0, q0);
diagrams.redraw; hold off;
