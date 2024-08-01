function [e_vec, R_0T_vec, p_0T_vec] = err_given_q12A(obj, qA_1, qA_2, qA_3, qB_3, qC_3, qD_3, qE_3, qF_3)
    e_vec = NaN([1 8]);
    R_0T_vec = NaN([3 3 8]);
    p_0T_vec = NaN([3 8]);

    % Given prismatic joints, find p_26 for each leg
    p_26_A = obj.kinA.P(:,3) + obj.kinA.P(:,4) + qA_3 * obj.kinA.H(:,3);
    p_26_B = obj.kinB.P(:,3) + obj.kinB.P(:,4) + qB_3 * obj.kinB.H(:,3);
    p_26_C = obj.kinC.P(:,3) + obj.kinC.P(:,4) + qC_3 * obj.kinC.H(:,3);
    p_26_D = obj.kinD.P(:,3) + obj.kinD.P(:,4) + qD_3 * obj.kinD.H(:,3);
    p_26_E = obj.kinE.P(:,3) + obj.kinE.P(:,4) + qE_3 * obj.kinE.H(:,3);
    p_26_F = obj.kinF.P(:,3) + obj.kinF.P(:,4) + qF_3 * obj.kinF.H(:,3);

    % 2D search over position of leg A
    R_02 = rot(obj.kinA.H(:,1), qA_1) * rot(obj.kinA.H(:,2), qA_2);
    p_06_A = obj.kinA.P(:,1) + R_02 * p_26_A;

    % Parameterize possible pose using constraint from leg B
    
    % Use Subproblem 3 to find alpha_B
    h_beta_B = normalize(p_06_A - obj.kinB.P(:,1));
    h_alpha_B = normalize(cross(p_26_B, h_beta_B )); % Careful
    [t_alpha_B, alpha_B_is_LS] = subproblem.sp_3(p_26_B, p_06_A - obj.kinB.P(:,1), h_alpha_B, norm(obj.kinB.P(:,end)-obj.kinA.P(:,end)));
    p_26_B_tilde = rot(h_alpha_B, t_alpha_B(1)) * p_26_B;
    
    % Use Subproblem 1 to find R1
    k2 = normalize(p_26_B_tilde + obj.kinB.P(:,1) - p_06_A);
    k1 = normalize(cross(k2, obj.kinA.P(:,end)-obj.kinB.P(:,end))); % Careful
    [theta1, theta1_is_LS] = subproblem.sp_1(obj.kinA.P(:,end)-obj.kinB.P(:,end), p_26_B_tilde + obj.kinB.P(:,1) - p_06_A, k1);
    R1 = rot(k1, theta1);
    
    % Find pose using 8th order subproblem and constraints from legs C, D
    [t_theta2, t_beta_B] = subproblem_8th_order.two_spheres(...
        R1*(obj.kinA.P(:,end)-obj.kinC.P(:,end)), obj.kinC.P(:,1)-p_06_A, ...
        R1*(obj.kinA.P(:,end)-obj.kinD.P(:,end)), obj.kinD.P(:,1)-p_06_A, ...
        k2, -h_beta_B, norm(p_26_C), norm(p_26_D));

    for i = 1:8
        theta_2 = t_theta2(i);
        beta_B = t_beta_B(i);
        % just use real solutions for now
        if imag(theta_2) ~= 0 || imag(beta_B) ~= 0
            continue;
        end

        R_06 = rot(h_beta_B, beta_B) * rot(k2, theta_2) * R1;
        p_0T = p_06_A + R_06 * obj.kinA.P(:,end);

        % Find error for lengths of legs E and F
        p_06_E = p_0T - R_06*obj.kinE.P(:,end);
        p_06_F = p_0T - R_06*obj.kinF.P(:,end);

        qE_3_i = norm(p_06_E - obj.kinE.P(:,1)) - norm(obj.kinE.P(:,3) + obj.kinE.P(:,4));
        qF_3_i = norm(p_06_F - obj.kinF.P(:,1)) - norm(obj.kinE.P(:,3) + obj.kinF.P(:,4));

        err_E = qE_3_i - qE_3; 
        err_F = qF_3_i - qF_3;

        R_0T_vec(:,:,i) = R_06;
        p_0T_vec(:,i) = p_0T;
        e_vec(:,i) = norm([err_E, err_F]);
    end

function e = normalize(v)
    e = v / norm(v);
end
end
