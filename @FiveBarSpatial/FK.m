function [R_0T_vec, p_0T_vec] = FK(obj, qA_active, qB_active, show_graph)
R_0T_vec = [];
p_0T_vec = [];


kinA = obj.kinA;
kinB = obj.kinB;

R_03_A = rot(kinA.H(:,1), qA_active(1))*rot(kinA.H(:,2), qA_active(2))*rot(kinA.H(:,3), qA_active(3));
R_03_B = rot(kinB.H(:,1), qB_active(1))*rot(kinB.H(:,2), qB_active(2))*rot(kinB.H(:,3), qB_active(3));

[q5A_vec, soln_num_vec] = search_1D( ...
@(q5A)(err_given_q5_A(q5A, kinA, kinB, R_03_A, R_03_B)), ...
-pi, pi, 2000, show_graph);

for i_q5A = 1:length(q5A_vec)
    [~, Q_partial_A, Q_partial_B, p_0T_vec_search] = err_given_q5_A(q5A_vec(i_q5A), kinA, kinB, R_03_A, R_03_B);
    q_partial_A = Q_partial_A(:, soln_num_vec(i_q5A));
    q_partial_B = Q_partial_B(:, soln_num_vec(i_q5A));

    R_05_A = R_03_A*rot(kinA.H(:,4), q_partial_A(1))*rot(kinA.H(:,5), q_partial_A(2));
    R_05_B = R_03_B*rot(kinB.H(:,4), q_partial_B(1))*rot(kinB.H(:,5), q_partial_B(2));
    % Find q6A using Subproblem 1
    [q6A, q6A_is_LS] = subproblem.sp_1(R_05_A'*R_05_B*kinB.H(:,6), kinB.R_6T'*kinB.H(:,6), -kinA.H(:,6));

    % Calculate R_0T
    R_0T_i = R_05_A*rot(kinA.H(:,6), q6A);

    R_0T_vec(:,:,i_q5A) = R_0T_i;
    p_0T_vec(:,i_q5A) = p_0T_vec_search(:, soln_num_vec(i_q5A));
end

end

function [e_vec, Q_partial_A, Q_partial_B, p_0T_vec] = err_given_q5_A(q5A, kinA, kinB, R_03_A, R_03_B)
    e_vec = NaN([1 4]);
    Q_partial_A = NaN([2 4]);
    Q_partial_B = NaN([2 4]);
    p_0T_vec = NaN([3 4]);
    i_soln = 1;

    R_45_A = rot(kinA.H(:,5), q5A);

    % Find q4A using Subproblem 3
    [t4A, t4A_is_LS] = subproblem.sp_3(...
        R_45_A*kinA.P(:,6),...
        R_03_A'*(kinB.P(:,1)+R_03_B*kinB.P(:,4))-kinA.P(:,4),...
        kinA.H(:,4),...
        norm(kinB.P(:,6)));
    if t4A_is_LS
        return;
    end

    for i_t4A = 1:length(t4A)
        q4A = t4A(i_t4A);
        % Calculate p_0T
        R_35_A = rot(kinA.H(:,4), q4A)*R_45_A;
        p_0T = R_03_A*(kinA.P(:,4)+R_35_A*kinA.P(:,6));

        % Find (q4B, q5B) using Subproblem 2
        [t4B, t5B, t45B_is_LS] = subproblem.sp_2(R_03_B'*(p_0T-kinB.P(:,1))-kinB.P(:,4), kinB.P(:,6), -kinB.H(:,4), kinB.H(:,5));
        if t45B_is_LS
            i_soln = i_soln + 2;
            continue;
        end

        for i_t45B = 1:length(t4B)
            q4B = t4B(i_t45B);
            q5B = t5B(i_t45B);
            
            R_05_A = R_03_A*R_35_A;
            R_05_B = R_03_B*rot(kinB.H(:,4), q4B)*rot(kinB.H(:,5), q5B);
            % Calculate error
            e_i = kinB.H(:,6)'*R_05_B'*R_05_A*kinA.H(:,6) - kinB.H(:,6)'*kinB.R_6T*kinB.H(:,6);
            e_vec(i_soln) = e_i;
            Q_partial_A(:,i_soln) = [q4A; q5A];
            Q_partial_B(:,i_soln) = [q4B; q5B];
            p_0T_vec(:,i_soln) = p_0T;
            i_soln = i_soln + 1;
        end
    end
end