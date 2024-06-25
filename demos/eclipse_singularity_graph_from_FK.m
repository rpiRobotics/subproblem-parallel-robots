rob = Eclipse;

% Constant p_0T
p_0T = [0;0;0];

% Vary R_0T
ex = [1;0;0];
ey = [0;1;0];

N_1 = 10;
N_2 = 10;
roll_vec  = linspace(-pi, pi, N_1);
pitch_vec = linspace(pi, pi, N_2);

for i_1 = 4%:N_1
for i_2 = 4%:N_2
    roll = roll_vec(i_1);
    pitch = pitch_vec(i_2);
    R_0T = rot(ey, pitch)*rot(ex, roll);

    % First IK
    [Q_A, Q_B, Q_C, is_LS_A, is_LS_B, is_LS_C] = rob.IK(R_0T, p_0T);
    q_A = Q_A(:,3);
    q_B = Q_B(:,3);
    q_C = Q_C(:,3);

    % Then FK
    [R_0T_vec, p_0T_vec, q_A_vec, q_B_vec, q_C_vec, slopes] = rob.FK(q_A(1:2), q_B(1:2), q_C(1:2))
    xline(q_A(3), 'r');
end
end