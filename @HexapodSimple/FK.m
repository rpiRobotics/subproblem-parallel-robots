function [R_0T_vec, p_0T_vec] = FK(obj, qA_3, qB_3, qC_3, qD_3, qE_3, qF_3, show_plot)
    R_0T_vec = [];
    p_0T_vec = [];

    % Given prismatic joints, find p_26 for each leg
    p_26_A = obj.kinA.P(:,3) + obj.kinA.P(:,4) + qA_3 * obj.kinA.H(:,3);
    p_26_B = obj.kinB.P(:,3) + obj.kinB.P(:,4) + qB_3 * obj.kinB.H(:,3);
    p_26_C = obj.kinC.P(:,3) + obj.kinC.P(:,4) + qC_3 * obj.kinC.H(:,3);
    p_26_D = obj.kinD.P(:,3) + obj.kinD.P(:,4) + qD_3 * obj.kinD.H(:,3);
    p_26_E = obj.kinE.P(:,3) + obj.kinE.P(:,4) + qE_3 * obj.kinE.H(:,3);
    p_26_F = obj.kinF.P(:,3) + obj.kinF.P(:,4) + qF_3 * obj.kinF.H(:,3);

    % Rotation axes for each of 3 circles formed by (A,B), (C,D), (E,F)
    h_beta_A = normalize(obj.kinB.P(:,1) - obj.kinA.P(:,1));
    h_beta_C = normalize(obj.kinD.P(:,1) - obj.kinC.P(:,1));
    h_beta_E = normalize(obj.kinF.P(:,1) - obj.kinE.P(:,1));

    % One point on circle for each of 3 circles using Subproblem 3
    % Careful! p_26_A and k_AB can't be parallel etc.
    h_alpha_A = normalize(cross(h_beta_A, p_26_A));
    h_alpha_C = normalize(cross(h_beta_C, p_26_C));
    h_alpha_E = normalize(cross(h_beta_E, p_26_E));

    [t_alpha_A] = subproblem.sp_3(p_26_A, obj.kinB.P(:,1) - obj.kinA.P(:,1), h_alpha_A, norm(p_26_B));
    [t_alpha_C] = subproblem.sp_3(p_26_C, obj.kinD.P(:,1) - obj.kinC.P(:,1), h_alpha_C, norm(p_26_D));
    [t_alpha_E] = subproblem.sp_3(p_26_E, obj.kinF.P(:,1) - obj.kinE.P(:,1), h_alpha_E, norm(p_26_F));

    % Subproblem 3 returns up to 2 results, but we only need to pick 1
    p_26_A_tilde = rot(h_alpha_A, t_alpha_A(1)) * p_26_A;
    p_26_C_tilde = rot(h_alpha_C, t_alpha_C(1)) * p_26_C;
    p_26_E_tilde = rot(h_alpha_E, t_alpha_E(1)) * p_26_E;

    % 1D search to find beta_A
    [beta_A_vec, soln_num_vec] = search_1D( ...
        @(x)err_given_beta_A(x, obj, h_beta_A, p_26_A_tilde, h_beta_C, p_26_C_tilde, h_beta_E, p_26_E_tilde), ...
        -pi, pi, 1000, show_plot);

    for i = 1:length(beta_A_vec)
        beta_A = beta_A_vec(i);
        soln_num = soln_num_vec(i);
        [~, p_06_A, p_06_C_vec, p_06_E_vec] = err_given_beta_A(beta_A, obj, h_beta_A, p_26_A_tilde, h_beta_C, p_26_C_tilde, h_beta_E, p_26_E_tilde);
        p_06_C = p_06_C_vec(:,soln_num);
        p_06_E = p_06_E_vec(:,soln_num);

        h_alpha = normalize(obj.kinA.P(:,end) - obj.kinC.P(:,end));
        h_gamma = normalize(p_06_C - p_06_A);
        h_beta = normalize(cross(h_alpha, h_gamma)); % Careful!

        % Find beta with Subproblem 1
        beta = subproblem.sp_1(obj.kinA.P(:,end) - obj.kinC.P(:,end), p_06_C - p_06_A, h_beta);
        R_beta = rot(h_beta, beta);

        % Find alpha with Subproblem 4
        % Keep only one solution
        t_alpha = subproblem.sp_4(R_beta'*h_gamma, obj.kinA.P(:,end) - obj.kinE.P(:,end), h_alpha, h_gamma'*(p_06_E - p_06_A));
        alpha = t_alpha(1);
        R_alpha = rot(h_alpha, alpha);

        % Find gamma with Subproblem 1
        gamma = subproblem.sp_1(R_beta*R_alpha*(obj.kinA.P(:,end) - obj.kinE.P(:,end)), p_06_E - p_06_A, h_gamma);
        R_gamma = rot(h_gamma, gamma);

        % Compute p_0T and R_0T
        R_0T = R_gamma*R_beta*R_alpha;
        p_0T = p_06_A + R_0T*obj.kinA.P(:,end);

        R_0T_vec(:,:,i) = R_0T;
        p_0T_vec(:,i) = p_0T;
    end
end

function [e_vec, p_06_A, p_06_C_vec, p_06_E_vec] = err_given_beta_A(beta_A, obj, h_beta_A, p_26_A_tilde, h_beta_C, p_26_C_tilde, h_beta_E, p_26_E_tilde)
    % Point on circle AB
    p_06_A = obj.kinA.P(:,1) + rot(h_beta_A, beta_A)*p_26_A_tilde;

    % Subproblem 3 to find beta_C
    [t_beta_C] = subproblem.sp_3(p_26_C_tilde, p_06_A-obj.kinC.P(:,1), h_beta_C, norm(obj.kinC.P(:,end) - obj.kinA.P(:,end)));
    % Duplicate if only 1 soln
    t_beta_C = [t_beta_C(1) t_beta_C(end)];

    p_06_C_1 = obj.kinC.P(:,1) + rot(h_beta_C, t_beta_C(1)) * p_26_C_tilde;
    p_06_C_2 = obj.kinC.P(:,1) + rot(h_beta_C, t_beta_C(2)) * p_26_C_tilde;

    % Subproblem 3 to find beta_E
    [t_beta_E] = subproblem.sp_3(p_26_E_tilde, p_06_A-obj.kinE.P(:,1), h_beta_E, norm(obj.kinE.P(:,end) - obj.kinA.P(:,end)));
    % Duplicate if only 1 soln
    t_beta_E = [t_beta_E(1) t_beta_E(end)];
    p_06_E_1 = obj.kinE.P(:,1) + rot(h_beta_E, t_beta_E(1)) * p_26_E_tilde;
    p_06_E_2 = obj.kinE.P(:,1) + rot(h_beta_E, t_beta_E(2)) * p_26_E_tilde;

    % Error in distance for p_06_C and p_06_E
    e_1 = norm(p_06_E_1 - p_06_C_1) - norm(obj.kinE.P(:,end) - obj.kinC.P(:,end));
    e_2 = norm(p_06_E_2 - p_06_C_1) - norm(obj.kinE.P(:,end) - obj.kinC.P(:,end));
    e_3 = norm(p_06_E_1 - p_06_C_2) - norm(obj.kinE.P(:,end) - obj.kinC.P(:,end));
    e_4 = norm(p_06_E_2 - p_06_C_2) - norm(obj.kinE.P(:,end) - obj.kinC.P(:,end));
    e_vec = [e_1; e_2; e_3; e_4];
    p_06_C_vec = [p_06_C_1 p_06_C_1 p_06_C_2 p_06_C_2];
    p_06_E_vec = [p_06_E_1 p_06_E_2 p_06_E_1 p_06_E_2];
end

function e = normalize(v)
    e = v / norm(v);
end