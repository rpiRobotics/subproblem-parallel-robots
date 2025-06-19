rob = DeltaBot;

% p_0T = rand_vec
p_0T = [0.4;0.1;-1];
% p_0T = [0;0;-1];

[Q_A, Q_B, Q_C, is_LS_A, is_LS_B, is_LS_C] = rob.IK(p_0T)
%%
i = 4
[R_A, T_A] = fwdkin(rob.kinA, Q_A(:,i))
[R_B, T_B] = fwdkin(rob.kinB, Q_B(:,i))
[R_C, T_C] = fwdkin(rob.kinC, Q_C(:,i))

%%
i = 4;
diagrams.setup; hold on
rob.plot(Q_A(:,i), Q_B(:,i), Q_C(:,i));
% diagrams.arrow([0;0;0], p_0T);

view(0, 30);

diagrams.redraw; hold off
%%
diagrams.save()

%% Test FK
i=4
[P_0T, is_LS_vec] = rob.FK(Q_A(1,i), Q_B(1,i), Q_C(1,i))

%%

% p_0T = P_0T(:,1); i = 1;
p_0T = P_0T(:,2); i = 4;


[Q_A, Q_B, Q_C, is_LS_A, is_LS_B, is_LS_C] = rob.IK(p_0T);
Q_A(:,i)
diagrams.setup; hold on
rob.plot(Q_A(:,i), Q_B(:,i), Q_C(:,i));

campos([0.2147  -11.3629    6.0174])
camtarget([ 0.2147         0   -0.5430])
camva(11.6435)

diagrams.redraw; hold off