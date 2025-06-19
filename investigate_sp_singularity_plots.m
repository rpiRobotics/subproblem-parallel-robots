% Subproblem 4 right at a singularity
h = [-1 0 0]';
p = [0 1 0]';
k = [0 0 1]';
d = 1;

[theta, is_LS] = subproblem.sp_4(h, p, k, d);
%%
search_1D(@(x)(search_err(x, h,p,k,d)), -pi, pi, 100, true);

function e = search_err(q1, h,p,k,d)
    e = norm(h'*rot(k,q1)*p - d)
end