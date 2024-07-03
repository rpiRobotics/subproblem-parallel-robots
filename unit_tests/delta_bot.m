rob = DeltaBot;

% p_0T = rand_vec
p_0T = [0.4;0.1;-1];

[Q_A, Q_B, Q_C, is_LS_A, is_LS_B, is_LS_C] = rob.IK(p_0T)

%%
i = 4;
diagrams.setup; hold on
rob.plot(Q_A(:,i), Q_B(:,i), Q_C(:,i));
diagrams.arrow([0;0;0], p_0T);

view(0, 30);

diagrams.redraw; hold off