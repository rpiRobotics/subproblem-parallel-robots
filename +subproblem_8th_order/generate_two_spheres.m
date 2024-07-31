% Generate coefficients for subproblem_8th_order.dual_quadratic

syms x1 x2 real
syms p1 p2 p3 p4 [3 1] real
syms k1 k2 [3 1] real
syms d1 d2 real

R1 = half_tan_rot(k1, x1);
R2 = half_tan_rot(k2, x2);

E1 = dot(p1, p1)/2 + dot(p2, p2)/2 - p1'*R1'*R2*p2 - d1^2/2
E2 = dot(p3, p3)/2 + dot(p4, p4)/2 - p3'*R1'*R2*p4 - d2^2/2

num1 = numden(E1);
num2 = numden(E2);

Pabc = fliplr(coeffs(num1, x2));
Qabc = fliplr(coeffs(num2, x2));

Pa = Pabc(1);
Pb = Pabc(2);
Pc = Pabc(3);
Qa = Qabc(1);
Qb = Qabc(2);
Qc = Qabc(3);

pa123 = fliplr(coeffs(Pa, x1));
pb123 = fliplr(coeffs(Pb, x1));
pc123 = fliplr(coeffs(Pc, x1));
qa123 = fliplr(coeffs(Qa, x1));
qb123 = fliplr(coeffs(Qb, x1));
qc123 = fliplr(coeffs(Qc, x1));


function R = half_tan_rot(k, x)
    s = 2*x/(x^2+1);
    c = (-x^2+1)/(x^2+1);
    R = k*k' + s*hat(k) - c*hat(k)^2;
end