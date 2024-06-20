function [Q, is_LS_vec, diag_sp_3, internal_1_vec, internal_2_vec] = IK(obj, T)
kin = obj.kin;

Q = [];
is_LS_vec = [];
internal_1_vec = [];
internal_2_vec = [];

T = T - kin.P(:,1);

[t2, t2_is_LS, diag_sp_3] = subproblem_comp.sp_3(kin.P(:,3), -kin.P(:,2), kin.H(:,2), norm(T));


for i_t2 = 1:length(t2)
    q2 = t2(i_t2);
    R_12 = rot(kin.H(:,2), q2);
    [q1, q1_is_LS, diag] = subproblem_comp.sp_1(kin.P(:,2)+R_12*kin.P(:,3), T, kin.H(:,1));
    

    internal_1_vec = [internal_1_vec diag.norm_k_X_p1];
    internal_2_vec = [internal_2_vec diag.norm_k_X_p1];
    Q = [Q [q1; q2]];
    is_LS_vec = [is_LS_vec t2_is_LS || q1_is_LS];
end


end