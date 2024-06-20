function [theta, is_LS, diag] = sp_1(p1, p2, k, theta_given)
arguments
    p1
    p2
    k
    theta_given = NaN;
end
% diagnostics
%   norm_k_X_p1 = 0 -> Internal singularity
%   norm_k_X_p2 = 0 -> Internal singularity

diag.norm_k_X_p1 = norm(cross(k, p1));
diag.norm_k_X_p1 = norm(cross(k, p2));

KxP = cross(k, p1);
A = [KxP -cross(k,KxP)];
pinv_A = A' / dot(KxP,KxP);

x = pinv_A*p2;

theta = atan2(x(1),x(2));

% Least squares solution if ||p_1|| != ||p_2|| or k'p_1 != k'p_2
is_LS = abs(norm(p1) - norm(p2)) > 1e-8 || abs(dot(k,p1) - dot(k,p2)) > 1e-8;

if ~isnan(theta_given)
    theta = theta_given;
end
end