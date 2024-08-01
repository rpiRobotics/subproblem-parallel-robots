function [R_0T_vec, p_0T_vec, qA_1_vec, qA_2_vec] = FK(obj, qA_3, qB_3, qC_3, qD_3, qE_3, qF_3, show_plot)

    [qA_1_vec, qA_2_vec, soln_num_vec] = search_2D(...
        @(qA_1, qA_2)obj.err_given_q12A(qA_1, qA_2, qA_3, qB_3, qC_3, qD_3, qE_3, qF_3),...
        -pi/2, pi/2, -pi, pi, 400, show_plot);
    
    N_solns = length(soln_num_vec);
    R_0T_vec = NaN(3, 3, N_solns);
    p_0T_vec = NaN(3, N_solns);

    for i = 1:length(soln_num_vec)
        qA_1_i = qA_1_vec(i);
        qA_2_i = qA_2_vec(i);
        [~, R_0T_vec_i, p_0T_vec_i] = obj.err_given_q12A(qA_1_i, qA_2_i, qA_3, qB_3, qC_3, qD_3, qE_3, qF_3);
        R_0T_vec(:,:,i) = R_0T_vec_i(:, :, soln_num_vec(i));
        p_0T_vec(:,i)   = p_0T_vec_i(:, soln_num_vec(i));
    end
end

