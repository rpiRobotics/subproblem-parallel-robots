function [theta, is_LS, diag] = sp_3(p1, p2, k, d)

[theta, is_LS, diag] = subproblem_comp.sp_4(p2, p1, k, 1/2 * (dot(p1,p1)+dot(p2,p2)-d^2));
end