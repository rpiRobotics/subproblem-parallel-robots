function [R_0T_vec, p_0T_vec, q_A_vec, q_B_vec, q_C_vec, slopes] = FK(obj, qA_active, qB_active, qC_active)
R_0T_vec = NaN(3,3,16);
p_0T_vec = NaN(3,1,16);
q_A_vec = NaN(6,2,16);
q_B_vec = NaN(6,2,16);
q_C_vec = NaN(6,2,16);

kinA = obj.kinA;
kinB = obj.kinB;
kinC = obj.kinC;

p_AB = kinA.R_6T'*kinA.P(:,end) - kinB.R_6T'*kinB.P(:,end);
p_BC = kinB.R_6T'*kinB.P(:,end) - kinC.R_6T'*kinC.P(:,end);
p_AC = kinA.R_6T'*kinA.P(:,end) - kinC.R_6T'*kinC.P(:,end);
k_AB = p_AB / norm(p_AB);
k_2 = cross(k_AB, p_BC); k_2 = k_2 / norm(k_2);
k_3 = cross(k_2, k_AB);

d_AB = norm(p_AB);
d_BC = norm(p_BC);
d_AC = norm(p_AC);

[q3_A_vec, soln_num_vec] = search_1D(@error_given_q3_A, -pi, pi, 1000, false);
slopes = subproblem_comp.get_1D_search_slope(@error_given_q3_A, q3_A_vec, soln_num_vec);

for i = 1:length(q3_A_vec)
    q3_A = q3_A_vec(i);
    soln_num = soln_num_vec(i);
    [~, p_06_A, p_06_B_vec, p_06_C_vec, q3_B_vec, q3_C_vec] = error_given_q3_A(q3_A);
    q3_B = q3_B_vec(soln_num);
    q3_C = q3_C_vec(soln_num);

    % Use Subproblem 2 to find theta_1 and theta_2
    [t2, t3, t_23_is_LS] = subproblem.sp_2(p_AB, p_06_B_vec(:,soln_num)-p_06_A,k_2, -k_3);
    theta_2 = t2(1); theta_3 = t3(1);
    R_2 = rot(k_2, theta_2);
    R_3 = rot(k_3, theta_3);

    % Use Subproblem 1 to find theta_3
    [t3, t3_is_LS] = subproblem.sp_1(p_AC, (R_3*R_2)'*(p_06_C_vec(:,soln_num)-p_06_A), k_AB);
    R_AB = rot(k_AB, t3);

    R_0T_i = R_3*R_2*R_AB;
    p_0T_i = p_06_A + R_0T_i * kinA.R_6T' * kinA.P(:,end);

    R_0T_vec(:,:,i) = R_0T_i;
    p_0T_vec(:,:,i) = p_0T_i;


    q_A_vec(:,:,i) = [qA_active(1:2) qA_active(1:2); q3_A q3_A; q_456(kinA, qA_active(1), q3_A, R_0T_i)];
    q_B_vec(:,:,i) = [qB_active(1:2) qB_active(1:2); q3_B q3_B; q_456(kinB, qB_active(1), q3_B, R_0T_i)];
    q_C_vec(:,:,i) = [qC_active(1:2) qC_active(1:2); q3_C q3_C; q_456(kinC, qC_active(1), q3_C, R_0T_i)];
end

function Q = q_456(kin, q_1, q_3, R_0T)
    Q = NaN(3,2);
    % Find (q_4, q_5) with Subproblem 2
    R_01 = rot(kin.H(:,1), q_1);
    R_23 = rot(kin.H(:,3), q_3);
    R_06 = R_0T * kin.R_6T';  
    R_36 = (R_01*R_23)'*R_06;
    [t5, t4, t45_is_LS] = subproblem.sp_2(kin.H(:,6), R_36*kin.H(:,6), kin.H(:,5), -kin.H(:,4));
    
    for i_45 = 1:length(t4)
        q4 = t4(i_45); q5 = t5(i_45);
        R_34 = rot(kin.H(:,4), q4);
        R_45 = rot(kin.H(:,5), q5);
        % Find q6 with Subproblem 1
        p = kin.H(:,5); % Not collinear with h_6;
        [q6, q6_is_LS] = subproblem.sp_1(p, (R_34*R_45)'*R_36*p, kin.H(:,6));
    
        Q(:,i_45) = [q4; q5; q6];
    end
end


function [e_vec, p_06_A, p_06_B_vec, p_06_C_vec, q3_B_vec, q3_C_vec] = error_given_q3_A(q3_A)
    e_vec = NaN(1,4);
    p_06_B_vec = NaN(3,16);
    p_06_C_vec = NaN(3,16);
    q3_B_vec = NaN(1,16);
    q3_C_vec = NaN(1,16);

    R_01_A = rot(kinA.H(:,1), qA_active(1));
    R_23_A = rot(kinA.H(:,3), q3_A);
    p_06_A = kinA.P(:,1) ...
        + R_01_A*(kinA.P(:,2) + kinA.P(:,3) + qA_active(2)*kinA.H(:,2))...
        + R_01_A*R_23_A*kinA.P(:,4);

    % Solve for q3_B using subproblem 3
    R_01_B = rot(kinB.H(:,1), qB_active(1));
    minus_p2_B = R_01_B'*(kinB.P(:,1)-p_06_A)+kinB.P(:,2)+kinB.P(:,3)+qB_active(2)*kinB.H(:,2);
    [t3_B, t3_B_is_LS] = subproblem.sp_3(kinB.P(:,4), -minus_p2_B, kinB.H(:,3), d_AB);
    if t3_B_is_LS
        return
    end

    % Solve for q3_C using subproblem 3
    R_01_C = rot(kinC.H(:,1), qC_active(1));
    minus_p2_C = R_01_C'*(kinC.P(:,1)-p_06_A)+kinC.P(:,2)+kinC.P(:,3)+qC_active(2)*kinC.H(:,2);
    [t3_C, t3_C_is_LS] = subproblem.sp_3(kinC.P(:,4), -minus_p2_C, kinC.H(:,3), d_AC);
    if t3_C_is_LS
        return
    end

    % Form error
    i_soln = 1;
    for i_t3_B = 1:length(t3_B)
    for i_t3_C = 1:length(t3_C)
        q3_B = t3_B(i_t3_B);
        q3_C = t3_C(i_t3_C);
        R_23_B = rot(kinB.H(:,3), q3_B);
        R_23_C = rot(kinC.H(:,3), q3_C);

        p_06_B = kinB.P(:,1) ...
        + R_01_B*(kinB.P(:,2) + kinB.P(:,3) + qB_active(2)*kinB.H(:,2))...
        + R_01_B*R_23_B*kinB.P(:,4);

        p_06_C = kinC.P(:,1) ...
        + R_01_C*(kinC.P(:,2) + kinC.P(:,3) + qC_active(2)*kinC.H(:,2))...
        + R_01_C*R_23_C*kinC.P(:,4);

        e_i = norm(p_06_C - p_06_B) - d_BC;
        e_vec(i_soln) = e_i;
        
        p_06_B_vec(:,i_soln) = p_06_B;
        p_06_C_vec(:,i_soln) = p_06_C;
        q3_B_vec(i_soln) = q3_B;
        q3_C_vec(i_soln) = q3_C;

        i_soln = i_soln + 1;
    end
    end
end
end