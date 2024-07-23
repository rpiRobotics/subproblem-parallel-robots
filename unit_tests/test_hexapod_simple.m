rob = HexapodSimple();

% [Q_A, Q_B, Q_C, Q_D, Q_E, Q_F] = rob.IK(eye(3), [0;0;2]);
[Q_A, Q_B, Q_C, Q_D, Q_E, Q_F] = rob.IK(rot([1;0;0], 0.4), [0;0;2]);
%%
i = 1;
diagrams.setup; hold on;
rob.plot(Q_A(:,i), Q_B(:,i), Q_C(:,i), Q_D(:,i), Q_E(:,i), Q_F(:,i));
diagrams.arrow([0;0;0], [0;0;2]);
diagrams.redraw; hold off;

%%
i = 1;
[R_0T_vec, p_0T_vec] = rob.FK(Q_A(3,i), Q_B(3,i), Q_C(3,i), Q_D(3,i), Q_E(3,i), Q_F(3,i), true);