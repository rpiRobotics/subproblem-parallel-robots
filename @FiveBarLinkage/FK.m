function [p_0T_vec, q_A_vec, q_B_vec] = FK(obj, q1_A, q1_B)
p_0T_vec = NaN(3,0);
q_A_vec = NaN(2,0);
q_B_vec = NaN(2,0);

kinA = obj.armA.kin;
kinB = obj.armB.kin;

R_01_A = rot(kinA.H(:,1), q1_A);
R_01_B = rot(kinB.H(:,1), q1_B);

minus_p2 = kinA.P(:,2) + R_01_A'*(kinA.P(:,1) - kinB.P(:,1) - R_01_B*kinB.P(:,2));
[t2_A, t_2A_is_LS] = subproblem.sp_3(kinA.P(:,3), -minus_p2, kinA.H(:,2), norm(kinB.P(:,3)));

for i_t2_A = 1:length(t2_A)
    q2_A = t2_A(i_t2_A);
    R_12_A = rot(kinA.H(:,2), q2_A);
    p_0T_A = kinA.P(:,1) + R_01_A*kinA.P(:,2) +  R_01_A*R_12_A*kinA.P(:,3);
    p_02_B = kinB.P(:,1) + R_01_B*kinB.P(:,2);

    [q2_B, q2_B_is_LS] = subproblem.sp_1(kinB.P(:,2), R_01_B'*(p_0T_A-p_02_B), kinB.H(:,2));

    p_0T_vec(:,i_t2_A) = p_0T_A;
    q_A_vec(:,i_t2_A) = [q1_A; q2_A];
    q_B_vec(:,i_t2_A) = [q1_B; q2_B];
end

end