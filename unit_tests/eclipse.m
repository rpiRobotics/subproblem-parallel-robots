rob = Eclipse;

R_0T = eye(3);
p_0T = [1.1;0.1;1.1];

[Q_A, Q_B, Q_C, is_LS_A, is_LS_B, is_LS_C] = rob.IK(R_0T, p_0T)

%%
i_soln = 2;
is_LS_A(i_soln)
q_B = Q_B(:,i_soln)

kin = rob.kinB;
[R_06_test, p_0T_test] = fwdkin(kin, q_B)
R_0T_test = R_06_test * kin.R_6T

%%
q_A = Q_A(:, i_soln);
q_B = Q_B(:, i_soln);
q_C = Q_C(:, i_soln);

diagrams.setup
hold on
rob.plot(q_A, q_B, q_C)
hold off
diagrams.redraw

%%
[R_06_inc, p_0T_inc] = fwdkin(rob.kinA, [0.1 0.2 0.3 0.4 0.5 0.6]')
R_0T_inc = R_06_inc * rob.kinA.R_6T

[Q_A, ~, ~, is_LS_A, ~, ~] = rob.IK(R_0T_inc, p_0T_inc)

%%
kin = rob.kinA;
k = cross(kin.H(:,2), kin.H(:,3));
p_06 = p_0T_inc - R_06_inc*kin.P(:,end);

R_23 = rot(kin.H(:,3), 0.3);
R_01 = rot(kin.H(:,1), 0.1);
LHS = k'*R_23'*kin.P(:,4)
RHS = k'*R_01'*(p_06 -kin.P(:,1)) - k'*(kin.P(:,2)+kin.P(:,3))

%% Animate a path

R_0T = eye(3);
p_0T_start = [-0.9;-.3;1.1];
p_0T_end = [1.1;0.3;1.1];

N = 30;
lambda = linspace(0,1,N);
P_0T = p_0T_end.*lambda + p_0T_start.*(1-lambda);

i_soln_A = 2;
i_soln_B = 4;
i_soln_C = 5;

for i = 1:N
    p_0T = P_0T(:,i);
    [Q_A, Q_B, Q_C, is_LS_A, is_LS_B, is_LS_C] = rob.IK(R_0T, p_0T);
    q_A = Q_A(:, i_soln_A);
    q_B = Q_B(:, i_soln_B);
    q_C = Q_C(:, i_soln_C);

    assert(is_LS_A(:,i_soln_A)==0)
    assert(is_LS_B(:,i_soln_B)==0)
    assert(is_LS_C(:,i_soln_C)==0)
    
    h_fig = diagrams.setup([5 5]); axis on; grid on
    view(0,30)
    hold on
    rob.plot(q_A, q_B, q_C)

    diagrams.arrow([0;0;0], p_0T);
    diagrams.line(p_0T_start, p_0T_end, lineStyle=":")
    diagrams.plane([0.1;0;0], [0;0;1],[1;0;0], 2.5);
    hold off

    diagrams.redraw();
    % exportgraphics(gca,"eclipse_IK.gif","Append",i ~= 1)
end

%% Test FK
rob = Eclipse;
[R_06_inc, ] = fwdkin(rob.kinA, [0.1 0.2 0.3 0.4 0.5 0.6]')
R_0T_inc = R_06_inc * rob.kinA.R_6T

[Q_A, Q_B, Q_C, is_LS_A, is_LS_B, is_LS_C] = rob.IK(R_0T_inc, p_0T_inc)

q_A = Q_A(:,6);
q_B = Q_B(:,6);
q_C = Q_C(:,6);
%%
[R_0T_vec, p_0T_vec, q_A_vec, q_B_vec, q_C_vec] = rob.FK(q_A(1:2), q_B(1:2), q_C(1:2))

%% Plot out all solutions
clear im
for i = 1:length(p_0T_vec)
    for j = 1:2
        h_fig = diagrams.setup(); %axis on; grid on
        % view(0,30)
        hold on
        rob.plot(q_A_vec(:,j,i), q_B_vec(:,j,i), q_C_vec(:,j,i))
        hold off
        camtarget([0 0 1])
        camva(12)
        campos([0.1626  -16.0599   10.7276])
        
        diagrams.redraw();

        drawnow
        frame = getframe(h_fig);
        im{2*(i-1)+j} = frame2im(frame);
    end
end
%%
filename = "eclipse_FK_2.gif";
for idx = 1:length(im)
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
        imwrite(A,map,filename,"gif","LoopCount",Inf,"DelayTime",0.5);
    else
        imwrite(A,map,filename,"gif","WriteMode","append","DelayTime",0.5);
    end
end
%%
[q3_A_vec(5) 0.3]
[R_0T_vec(:,:,5) R_0T_inc]
[p_0T_vec(:,:,5) p_0T_inc]

%% Verify eqn
kinA = rob.kinA;
kinB = rob.kinB;
R_01_A = rot(kinA.H(:,1), q_A(1));
R_23_A = rot(kinA.H(:,3), q_A(3));
R_01_B = rot(kinB.H(:,1), q_B(1));
R_23_B = rot(kinB.H(:,3), q_B(3));

p_06_A = kinA.P(:,1) ...
    + R_01_A*(kinA.P(:,2) + kinA.P(:,3) + q_A(2)*kinA.H(:,2))...
    + R_01_A*R_23_A*kinA.P(:,4);

p_06_B = kinB.P(:,1) ...
    + R_01_B*(kinB.P(:,2) + kinB.P(:,3) + q_B(2)*kinB.H(:,2))...
    + R_01_B*R_23_B*kinB.P(:,4);

R_0T_inc*(kinA.R_6T'*kinA.P(:,end) - kinB.R_6T'*kinB.P(:,end))
p_06_B-p_06_A
