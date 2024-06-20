function [theta, is_LS, diag] = sp_4(h, p, k, d, theta_given)
arguments
    h
    p
    k
    d
    theta_given = NaN;
end

% diagnostics
%   norm_x = 1 -> boundary singularity
%   norm_k_X_p = 0 -> internal singularity, theta arbitrary
%   norm_k_X_h = 0 -> internal singularity, theta arbitrary

diag.norm_k_X_p = norm(cross(k,p));
diag.norm_k_X_h = norm(cross(k,h));

A_11 = cross(k,p);
A_1 = [A_11 -cross(k,A_11)];
A = h'*A_1;

b = d - h'*k*(k'*p);

norm_A_2 = dot(A,A);

x_min = A'/norm_A_2 * b;
diag.norm_x = norm(x_min);

if diag.norm_x < 1
    xi = sqrt(norm_A_2-b^2) / norm_A_2;
    x_N_prime = [A(2); -A(1)];

    sc_1 = x_min + xi*x_N_prime;
    sc_2 = x_min - xi*x_N_prime;

    theta = [atan2(sc_1(1), sc_1(2)) atan2(sc_2(1), sc_2(2))];
    is_LS = false;
else
    theta = atan2(x_min(1), x_min(2));
    is_LS = true;
end

if ~isnan(theta_given)
    theta = theta_given;
end

end