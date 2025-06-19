rob = Eclipse;
%%
% q_A = [0; 0];
% q_B = [deg2rad(120); 0.5];
% q_C = [deg2rad(240); 2.4298];

% q_A = [0; 0];
% q_B = [deg2rad(120); 0.2];
% q_C = [deg2rad(297.92); 0.4];

% q_A = [0; 0];
% q_B = [deg2rad(120); 0.2];
% q_C = [deg2rad(240); 0.31];

q_A = [0; 0];
q_B = [deg2rad(30); 0.2635];
q_C = [deg2rad(60); 0.2];

% q_A = [0; 0];
% q_B = [deg2rad(120); 0.5];
% q_C = [deg2rad(240); 3.21];
% [R_0T_vec, p_0T_vec, q_A_vec, q_B_vec, q_C_vec, slopes] = rob.FK(q_A(1:2), q_B(1:2), q_C(1:2))

%%
[R_0T_vec, p_0T_vec, q_A_vec, q_B_vec, q_C_vec, slopes] = rob.FK(q_C(1:2), q_A(1:2), q_B(1:2))
%%
% Make a graph of error zeros vs some paramter
N = 100;
q_vec = linspace(-pi, pi, N);
q3_A_MAT = NaN(N, 16);
for i = 1:N
    [R_0T_vec, p_0T_vec, q_A_vec, q_B_vec, q_C_vec] = rob.FK(q_A, [q_vec(i); q_B(2)], q_C);
    q3_A_i = q_A_vec(3,1,:);
    q3_A_i = q3_A_i(:);
    q3_A_MAT(i,1:length(q3_A_i)) = q3_A_i;
    i
end
%%
plot(q_vec, q3_A_MAT, '.k')
%%
for i = 3
    for j = 1
        h_fig = diagrams.setup(); %axis on; grid on
        annotation('rectangle',[0 0 1 1]);
        hold on
        rob.plot(q_A_vec(:,j,i), q_B_vec(:,j,i), q_C_vec(:,j,i))
        hold off
        % camtarget([-0.5 0 1.6])
        % camva(10)
        % campos([0  -15   10])
        diagrams.redraw();
        % exportgraphics(gca, "eclipse_FK_singularity.gif", "Append", i ~= 7)
    end
end

%% Check with Jacobian

% Unstable singularity if JT_p * JC_p_tilde =/= 0

% JT_p maps passive joints velocities to end effector velocity
% JC_p_tilde is the span of possible constrained passive joint velocities
% JC_p_tilde spans the null space of JC_p
% JC_p forms the constraint equation for the passive joint velocities

q_A = q_A_vec(:,j,i);
q_B = q_B_vec(:,j,i);
q_C = q_C_vec(:,j,i);

J_A = robotjacobian(rob.kinA, q_A);
J_B = robotjacobian(rob.kinA, q_B);
J_C = robotjacobian(rob.kinA, q_C);

J_A_p = J_A(:,3:6);
J_B_p = J_B(:,3:6);
J_C_p = J_C(:,3:6);
z64 = zeros(6,4);
z66 = zeros(6,6);

JT_p = [J_A_p J_B_p J_C_p];

% Constraint 1: vA = vB
% Constraint 2: vA = vC
JC_p = [J_A_p -J_B_p    z64
        J_A_p    z64 -J_C_p];
svd(JC_p);
JC_p_tilde = null(JC_p, 1e-1);
JC_p*JC_p_tilde

JT_p * JC_p_tilde
