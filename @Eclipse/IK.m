function [Q_A, Q_B, Q_C, is_LS_A, is_LS_B, is_LS_C] = IK(obj, R_0T, p_0T)
    [Q_A, is_LS_A] = IK_single_arm(obj.kinA, R_0T, p_0T);
    [Q_B, is_LS_B] = IK_single_arm(obj.kinB, R_0T, p_0T);
    [Q_C, is_LS_C] = IK_single_arm(obj.kinC, R_0T, p_0T);
end

function [Q, is_LS_vec] = IK_single_arm(kin, R_0T, p_0T)
Q = [];
is_LS_vec = [];

R_06 = R_0T * kin.R_6T';
p_06 = p_0T - R_06*kin.P(:,end);

k = cross(kin.H(:,2), kin.H(:,3));

% Solve for q_1 using Subproblem 4
[t1, t1_is_LS] = subproblem.sp_4(kin.H(:,3), p_06 - kin.P(:,1), -kin.H(:,1), ...
    dot(kin.H(:,3), kin.P(:,2)+kin.P(:,3)+kin.P(:,4)));
for i_t1 = 1:length(t1)
    q1 = t1(i_t1);
    R_01 = rot(kin.H(:,1), q1);
    % Solve for q_3 using Subproblem 4
    [t3, t3_is_LS] = subproblem.sp_4(k, kin.P(:,4), kin.H(:,3), ...
        k'*R_01'*(p_06-kin.P(:,1))-k'*(kin.P(:,2)+kin.P(:,3)));
    % Project to find q_2
    for i_t3 = 1:length(t3)
        q3 = t3(i_t3);
        R_23 = rot(kin.H(:,3), q3);
        q2 = dot(kin.H(:,2), ...
            p_06-kin.P(:,1)-kin.P(:,2)-kin.P(:,3)-R_23*kin.P(:,4));

        % Find (q_4, q_5) with Subproblem 2
        R_36 = (R_01*R_23)'*R_06;
        [t5, t4, t45_is_LS] = subproblem.sp_2(kin.H(:,6), R_36*kin.H(:,6), kin.H(:,5), -kin.H(:,4));
        for i_45 = 1:length(t4)
            q4 = t4(i_45); q5 = t5(i_45);
            R_34 = rot(kin.H(:,4), q4);
            R_45 = rot(kin.H(:,5), q5);
            % Find q6 with Subproblem 1
            p = kin.H(:,5); % Not collinear with h_6;
            [q6, q6_is_LS] = subproblem.sp_1(p, (R_34*R_45)'*R_36*p, kin.H(:,6));

            Q = [Q [q1;q2;q3;q4;q5;q6]];
            is_LS_vec = [is_LS_vec t1_is_LS||t3_is_LS||t45_is_LS||q6_is_LS];

        end


    end
end

end