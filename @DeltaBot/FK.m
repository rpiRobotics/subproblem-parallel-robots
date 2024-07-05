function [P_0T, is_LS_vec] = FK(obj, q1_A, q1_B, q1_C)

P_0T = [];
is_LS_vec = [];

kinA = obj.kinA;
kinB = obj.kinB;
kinC = obj.kinC;

p_02_Ap = kinA.P(:,1) + rot(kinA.H(:,1), q1_A)*kinA.P(:,2) + kinA.P(:,6);
p_02_Bp = kinB.P(:,1) + rot(kinB.H(:,1), q1_B)*kinB.P(:,2) + kinB.P(:,6);
p_02_Cp = kinC.P(:,1) + rot(kinC.H(:,1), q1_C)*kinC.P(:,2) + kinC.P(:,6);

h_beta = normalize(p_02_Bp - p_02_Ap);
h_alpha = normalize(cross([0;0;1], h_beta)); % Careful! Make sure h_beta is not parallel to [0;0;1]

p_34p = norm(kinA.P(:,4)) * h_beta;

% Find alpha with Subproblem 3.
[alpha_vec, alpha_is_LS] = subproblem.sp_3(p_34p, p_02_Bp - p_02_Ap, h_alpha, norm(kinB.P(:,4)));

% Keep only one value of alpha
alpha = alpha_vec(1);

% find beta with Subproblem 3
[beta_vec, beta_is_LS] = subproblem.sp_3(rot(h_alpha, alpha)*p_34p, p_02_Cp - p_02_Ap, h_beta, norm(kinC.P(:,4)));

for i_beta = 1:length(beta_vec)
    beta = beta_vec(i_beta);
    p_0T = p_02_Ap + rot(h_beta, beta)*rot(h_alpha, alpha)*p_34p;
    P_0T = [P_0T p_0T];
    % is_LS_vec = [is_LS_vec alpha_is_LS||beta_is_LS];
    is_LS_vec = [is_LS_vec [alpha_is_LS;beta_is_LS]];
end
end


function n = normalize(p)
    n = p / norm(p);
end