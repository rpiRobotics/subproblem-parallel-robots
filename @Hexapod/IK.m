function [Q_A, Q_B, Q_C, Q_D, Q_E, Q_F] = IK(obj, R_0T, p_0T)
    [Q_A] = IK_single_arm(obj.kinA, R_0T, p_0T);
    [Q_B] = IK_single_arm(obj.kinB, R_0T, p_0T);
    [Q_C] = IK_single_arm(obj.kinC, R_0T, p_0T);
    [Q_D] = IK_single_arm(obj.kinD, R_0T, p_0T);
    [Q_E] = IK_single_arm(obj.kinE, R_0T, p_0T);
    [Q_F] = IK_single_arm(obj.kinF, R_0T, p_0T);
end

function [Q, is_LS] = IK_single_arm(kin, R_0T, p_0T)
    Q = [];
    is_LS = [];
    R_06 = R_0T;
    p_06 = p_0T - R_06*kin.P(:,end) - kin.P(:,1);

    q3 = norm(p_06)-norm(kin.P(:,3)+kin.P(:,4));

    [t1, t2, t12_is_LS] = subproblem.sp_2(p_06, kin.P(:,3)+kin.P(:,4), -kin.H(:,1), kin.H(:,2));
    for i_12 = 1:length(t1)
        q1 = t1(i_12);
        q2 = t2(i_12);
        R_12 = rot(kin.H(:,1), q1)*rot(kin.H(:,2), q2);
        R_36 = R_12'*R_06;

        [t4, t5, t45_is_LS] = subproblem.sp_2(R_36*kin.H(:,6), kin.H(:,6), -kin.H(:,4), kin.H(:,5));
        for i_45 = 1:length(t4)
            q4 = t4(i_45);
            q5 = t5(i_45);
            R_05 = R_12*rot(kin.H(:,4), q4)*rot(kin.H(:,5), q5);
            [q6, q6_is_LS] = subproblem.sp_1(kin.H(:,5), R_05'*R_06*kin.H(:,5), kin.H(:,6));

            Q = [Q [q1;q2;q3;q4;q5;q6]];
            is_LS = [is_LS t12_is_LS||t45_is_LS||q6_is_LS];
        end
    end
end