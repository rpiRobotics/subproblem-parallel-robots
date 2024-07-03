function [Q_A, Q_B, Q_C, is_LS_A, is_LS_B, is_LS_C] = IK(obj, p_0T)
    [Q_A, is_LS_A] = IK_single_arm(obj.kinA, p_0T);
    [Q_B, is_LS_B] = IK_single_arm(obj.kinB, p_0T);
    [Q_C, is_LS_C] = IK_single_arm(obj.kinC, p_0T);
end

function [Q, is_LS_vec] = IK_single_arm(kin, p_0T)
Q = [];
is_LS_vec = [];

p_04 = p_0T - kin.R_5T'*kin.P(:,end);

[t1, t1_is_LS] = subproblem.sp_3(p_04 - kin.P(:,1), kin.P(:,2), -kin.H(:,1), norm(kin.P(:,4)));

for i_t1 = 1:length(t1)
    q1 = t1(i_t1);
    R_01 = rot(kin.H(:,1), q1);

    [t2, t3, t23_is_LS] = subproblem.sp_2(R_01'*(p_04 - kin.P(:,1)) - kin.P(:,2), kin.P(:,4), -kin.H(:,2), kin.H(:,3));

    for i_t23 = 1:length(t2)
        q2 = t2(i_t23);
        q3 = t3(i_t23);
        R_12 = rot(kin.H(:,2), q2);
        R_23 = rot(kin.H(:,3), q3);
        R_03 = R_01*R_12*R_23;

        % q4 should always equal -q3
        [q4, q4_is_LS] = subproblem.sp_1(kin.H(:,5), R_03'*kin.R_5T'*kin.H(:,5), kin.H(:,4));
        [q5, q5_is_LS] = subproblem.sp_1(kin.H(:,4), kin.R_5T*R_03  *kin.H(:,4), -kin.H(:,5));

        Q = [Q [q1;q2;q3;q4;q5]];
        is_LS_vec = [is_LS_vec t1_is_LS||t23_is_LS||q4_is_LS||q5_is_LS];
        % is_LS_vec = [is_LS_vec [t1_is_LS; t23_is_LS; q4_is_LS; q5_is_LS]];
    end
end

end