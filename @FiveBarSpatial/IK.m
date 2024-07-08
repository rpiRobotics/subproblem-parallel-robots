function [Q_A, Q_B, is_LS_A, is_LS_B] = IK(obj, R_0T, p_0T)
    [Q_A, is_LS_A] = IK_single_arm(obj.kinA, R_0T, p_0T);
    [Q_B, is_LS_B] = IK_single_arm(obj.kinB, R_0T, p_0T);
end

function [Q, is_LS] = IK_single_arm(kin, R_0T, p_0T)
    Q = [];
    is_LS = [];
    R_06 = R_0T*kin.R_6T';
    p_06 = p_0T - kin.P(:,1);

    % Solve for q6 with Subproblem 3
    [t6, t6_is_LS] = subproblem.sp_3(kin.P(:,6), R_06'*p_06, -kin.H(:,6), norm(kin.P(:,4)));
    
    for i_t6 = 1:length(t6)
        q6 = t6(i_t6);
        R_56 = rot(kin.H(:,6), q6);

         % Solve for (q4, q5) with Subproblem 2
        [t4, t5, t45_is_LS] = subproblem.sp_2(kin.P(:,4), R_56*R_06'*p_06 - kin.P(:,6), -kin.H(:,4), kin.H(:,5));

        for i_45 = 1:length(t4)
            q4 = t4(i_45);
            q5 = t5(i_45);
            R_46 = rot(kin.H(:,4), q4)*rot(kin.H(:,5), q5)*R_56;

            % Solve for spherical shoulder using Subproblems 4 and 1
            R_03  = R_06*R_46';

            % Solve for (q1, q2) with Subproblem 2
            [t1, t2, t12_is_LS] = subproblem.sp_2(R_03*kin.H(:,3), kin.H(:,3), -kin.H(:,1), kin.H(:,2));
            for i_12 = 1:length(t1)
                q1 = t1(i_12);
                q2 = t2(i_12);
                R_12 = rot(kin.H(:,1), q1)*rot(kin.H(:,2), q2);

                % Solve for q3 with Subproblem 1
                [q3, q3_is_LS] = subproblem.sp_1(kin.H(:,2), R_12'*R_03*kin.H(:,2), kin.H(:,3));

                Q = [Q [q1;q2;q3;q4;q5;q6]];
                is_LS = [is_LS t6_is_LS||t45_is_LS||t12_is_LS||q3_is_LS];
            end
        end
    end
end