rob = FiveBarLinkage;

p_0T = [0.5; 0.5; 0];

[Q_A, Q_B, is_LS_A, is_LS_B] =  rob.IK(p_0T);

q_A = Q_A(:,1);
q_B = Q_B(:,2);

diagrams.setup;
hold on
view(2)

rob.plot(q_A, q_B);

diagrams.arrow([0;0;0], p_0T);

hold off
diagrams.redraw;

%%
[p_0T_vec, q_A_vec, q_B_vec] = rob.FK(q_A(1), q_B(1))