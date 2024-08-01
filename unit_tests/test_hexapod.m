rob = Hexapod();
%% Each arm in zero pose
q0 = zeros(6,1);

diagrams.setup; hold on;
rob.plot(q0, q0, q0, q0, q0, q0);
diagrams.redraw; hold off;
%% Whole-arm IK
[Q_A, Q_B, Q_C, Q_D, Q_E, Q_F] = rob.IK(eye(3), [0;0;2]);

%%
R_0T = rot([1;0;0], 0.1)*rot([0;1;0], 0.2)*rot([0;0;1], 0.3);
p_0T = [0;0;2];

[Q_A, Q_B, Q_C, Q_D, Q_E, Q_F] = rob.IK(R_0T, p_0T);

%%
i = 1;
diagrams.setup; hold on;
rob.plot(Q_A(:,i), Q_B(:,i), Q_C(:,i), Q_D(:,i), Q_E(:,i), Q_F(:,i));
diagrams.arrow([0;0;0], [0;0;2]);
diagrams.redraw; hold off;

%% Single arm FK
[R_t, T_t] = fwdkin(rob.kinC, Q_C(:,i))

%% Test error function
i = 1;
[e_vec, R_0T_vec, p_0T_vec] = rob.err_given_q12A(Q_A(1,i), Q_A(2,i), Q_A(3,i), Q_B(3,i), Q_C(3,i), Q_D(3,i), Q_E(3,i), Q_F(3,i))

%% Test 2D search FK

[R_0T_vec, p_0T_vec, qA_1_vec, qA_2_vec] = rob.FK(Q_A(3,i), Q_B(3,i), Q_C(3,i), Q_D(3,i), Q_E(3,i), Q_F(3,i), true)


%% Verify solutions

i_IK = 7;

[R_0T_vec(:,:,i_IK), p_0T_vec(:,i_IK)]

[Q_A_t, Q_B_t, Q_C_t, Q_D_t, Q_E_t, Q_F_t] = rob.IK(R_0T_vec(:,:,i_IK), p_0T_vec(:,i_IK));

eA = abs(Q_A_t(3) -  Q_A(3,i));
eB = abs(Q_B_t(3) -  Q_B(3,i));
eC = abs(Q_C_t(3) -  Q_C(3,i));
eD = abs(Q_D_t(3) -  Q_D(3,i));
eE = abs(Q_E_t(3) -  Q_E(3,i));
eF = abs(Q_F_t(3) -  Q_F(3,i));

[eA eB eC eD eE eF]

%%
i_IK = 1;
[e_vec] = rob.err_given_q12A(qA_1_vec(i_IK), qA_2_vec(i_IK), Q_A(3,i), Q_B(3,i), Q_C(3,i), Q_D(3,i), Q_E(3,i), Q_F(3,i))