% Solve the following system of polynomials
% Pa(x1) x2^2 + Pb(x1) x2 + Pc(x1) = 0
% Qa(x1) x2^2 + Qb(x1) x2 + Qc(x1) = 0
%
% Each polynomial P and Q is quadratic
% Represent polynomials by vectors of coefficients of x1
% P = [p1 p2 p3] where p1 x1^2 + p2 x1 + p3 = 0

function [x1, x2] = dual_quadratic(Pa, Pb, Pc, Qa, Qb, Qc)

pa1 = Pa(1); pa2 = Pa(2); pa3 = Pa(3);
pb1 = Pb(1); pb2 = Pb(2); pb3 = Pb(3);
pc1 = Pc(1); pc2 = Pc(2); pc3 = Pc(3);
qa1 = Qa(1); qa2 = Qa(2); qa3 = Qa(3);
qb1 = Qb(1); qb2 = Qb(2); qb3 = Qb(3);
qc1 = Qc(1); qc2 = Qc(2); qc3 = Qc(3);


% These polynomial coefficients were generated using generate_dual_quadratic.m
C = NaN([9 1]);
C(9) = pa3^2*qc3^2 - pa3*pb3*qb3*qc3 - 2*pa3*pc3*qa3*qc3 + pa3*pc3*qb3^2 + pb3^2*qa3*qc3 - pb3*pc3*qa3*qb3 + pc3^2*qa3^2;
C(8) = 2*qc2*pa3^2*qc3 - qc2*pa3*pb3*qb3 - qb2*pa3*pb3*qc3 - 2*qc2*pa3*pc3*qa3 + 2*qb2*pa3*pc3*qb3 - 2*qa2*pa3*pc3*qc3 - 2*pc2*pa3*qa3*qc3 + pc2*pa3*qb3^2 - pb2*pa3*qb3*qc3 + 2*pa2*pa3*qc3^2 + qc2*pb3^2*qa3 + qa2*pb3^2*qc3 - qb2*pb3*pc3*qa3 - qa2*pb3*pc3*qb3 - pc2*pb3*qa3*qb3 + 2*pb2*pb3*qa3*qc3 - pa2*pb3*qb3*qc3 + 2*qa2*pc3^2*qa3 + 2*pc2*pc3*qa3^2 - pb2*pc3*qa3*qb3 - 2*pa2*pc3*qa3*qc3 + pa2*pc3*qb3^2;
C(7) = pa2^2*qc3^2 + 4*pa2*pa3*qc2*qc3 - pa2*pb2*qb3*qc3 - pa2*pb3*qb2*qc3 - pa2*pb3*qb3*qc2 - 2*pa2*pc2*qa3*qc3 + pa2*pc2*qb3^2 - 2*pa2*pc3*qa2*qc3 - 2*pa2*pc3*qa3*qc2 + 2*pa2*pc3*qb2*qb3 + pa3^2*qc2^2 + 2*qc1*pa3^2*qc3 - pa3*pb2*qb2*qc3 - pa3*pb2*qb3*qc2 - pa3*pb3*qb2*qc2 - qc1*pa3*pb3*qb3 - qb1*pa3*pb3*qc3 - 2*pa3*pc2*qa2*qc3 - 2*pa3*pc2*qa3*qc2 + 2*pa3*pc2*qb2*qb3 - 2*pa3*pc3*qa2*qc2 - 2*qc1*pa3*pc3*qa3 + pa3*pc3*qb2^2 + 2*qb1*pa3*pc3*qb3 - 2*qa1*pa3*pc3*qc3 - 2*pc1*pa3*qa3*qc3 + pc1*pa3*qb3^2 - pb1*pa3*qb3*qc3 + 2*pa1*pa3*qc3^2 + pb2^2*qa3*qc3 + 2*pb2*pb3*qa2*qc3 + 2*pb2*pb3*qa3*qc2 - pb2*pc2*qa3*qb3 - pb2*pc3*qa2*qb3 - pb2*pc3*qa3*qb2 + pb3^2*qa2*qc2 + qc1*pb3^2*qa3 + qa1*pb3^2*qc3 - pb3*pc2*qa2*qb3 - pb3*pc2*qa3*qb2 - pb3*pc3*qa2*qb2 - qb1*pb3*pc3*qa3 - qa1*pb3*pc3*qb3 - pc1*pb3*qa3*qb3 + 2*pb1*pb3*qa3*qc3 - pa1*pb3*qb3*qc3 + pc2^2*qa3^2 + 4*pc2*pc3*qa2*qa3 + pc3^2*qa2^2 + 2*qa1*pc3^2*qa3 + 2*pc1*pc3*qa3^2 - pb1*pc3*qa3*qb3 - 2*pa1*pc3*qa3*qc3 + pa1*pc3*qb3^2;
C(6) = 2*pa2^2*qc2*qc3 + 2*pa2*pa3*qc2^2 + 4*qc1*pa2*pa3*qc3 - pa2*pb2*qb2*qc3 - pa2*pb2*qb3*qc2 - pa2*pb3*qb2*qc2 - qc1*pa2*pb3*qb3 - qb1*pa2*pb3*qc3 - 2*pa2*pc2*qa2*qc3 - 2*pa2*pc2*qa3*qc2 + 2*pa2*pc2*qb2*qb3 - 2*pa2*pc3*qa2*qc2 - 2*qc1*pa2*pc3*qa3 + pa2*pc3*qb2^2 + 2*qb1*pa2*pc3*qb3 - 2*qa1*pa2*pc3*qc3 - 2*pc1*pa2*qa3*qc3 + pc1*pa2*qb3^2 - pb1*pa2*qb3*qc3 + 2*pa1*pa2*qc3^2 + 2*qc1*pa3^2*qc2 - pa3*pb2*qb2*qc2 - qc1*pa3*pb2*qb3 - qb1*pa3*pb2*qc3 - qc1*pa3*pb3*qb2 - qb1*pa3*pb3*qc2 - 2*pa3*pc2*qa2*qc2 - 2*qc1*pa3*pc2*qa3 + pa3*pc2*qb2^2 + 2*qb1*pa3*pc2*qb3 - 2*qa1*pa3*pc2*qc3 - 2*qc1*pa3*pc3*qa2 + 2*qb1*pa3*pc3*qb2 - 2*qa1*pa3*pc3*qc2 - 2*pc1*pa3*qa2*qc3 - 2*pc1*pa3*qa3*qc2 + 2*pc1*pa3*qb2*qb3 - pb1*pa3*qb2*qc3 - pb1*pa3*qb3*qc2 + 4*pa1*pa3*qc2*qc3 + pb2^2*qa2*qc3 + pb2^2*qa3*qc2 + 2*pb2*pb3*qa2*qc2 + 2*qc1*pb2*pb3*qa3 + 2*qa1*pb2*pb3*qc3 - pb2*pc2*qa2*qb3 - pb2*pc2*qa3*qb2 - pb2*pc3*qa2*qb2 - qb1*pb2*pc3*qa3 - qa1*pb2*pc3*qb3 - pc1*pb2*qa3*qb3 + 2*pb1*pb2*qa3*qc3 - pa1*pb2*qb3*qc3 + qc1*pb3^2*qa2 + qa1*pb3^2*qc2 - pb3*pc2*qa2*qb2 - qb1*pb3*pc2*qa3 - qa1*pb3*pc2*qb3 - qb1*pb3*pc3*qa2 - qa1*pb3*pc3*qb2 - pc1*pb3*qa2*qb3 + 2*pb1*pb3*qa2*qc3 - pc1*pb3*qa3*qb2 + 2*pb1*pb3*qa3*qc2 - pa1*pb3*qb2*qc3 - pa1*pb3*qb3*qc2 + 2*pc2^2*qa2*qa3 + 2*pc2*pc3*qa2^2 + 4*qa1*pc2*pc3*qa3 + 2*pc1*pc2*qa3^2 - pb1*pc2*qa3*qb3 - 2*pa1*pc2*qa3*qc3 + pa1*pc2*qb3^2 + 2*qa1*pc3^2*qa2 + 4*pc1*pc3*qa2*qa3 - pb1*pc3*qa2*qb3 - 2*pa1*pc3*qa2*qc3 - pb1*pc3*qa3*qb2 - 2*pa1*pc3*qa3*qc2 + 2*pa1*pc3*qb2*qb3;
C(5) = pa1^2*qc3^2 + 4*pa1*pa2*qc2*qc3 + 4*pa1*pa3*qc1*qc3 + 2*pa1*pa3*qc2^2 - pa1*pb1*qb3*qc3 - pa1*pb2*qb2*qc3 - pa1*pb2*qb3*qc2 - pa1*pb3*qb1*qc3 - pa1*pb3*qb2*qc2 - pa1*pb3*qb3*qc1 - 2*pa1*pc1*qa3*qc3 + pa1*pc1*qb3^2 - 2*pa1*pc2*qa2*qc3 - 2*pa1*pc2*qa3*qc2 + 2*pa1*pc2*qb2*qb3 - 2*pa1*pc3*qa1*qc3 - 2*pa1*pc3*qa2*qc2 - 2*pa1*pc3*qa3*qc1 + 2*pa1*pc3*qb1*qb3 + pa1*pc3*qb2^2 + 2*pa2^2*qc1*qc3 + pa2^2*qc2^2 + 4*pa2*pa3*qc1*qc2 - pa2*pb1*qb2*qc3 - pa2*pb1*qb3*qc2 - pa2*pb2*qb1*qc3 - pa2*pb2*qb2*qc2 - pa2*pb2*qb3*qc1 - pa2*pb3*qb1*qc2 - pa2*pb3*qb2*qc1 - 2*pa2*pc1*qa2*qc3 - 2*pa2*pc1*qa3*qc2 + 2*pa2*pc1*qb2*qb3 - 2*pa2*pc2*qa1*qc3 - 2*pa2*pc2*qa2*qc2 - 2*pa2*pc2*qa3*qc1 + 2*pa2*pc2*qb1*qb3 + pa2*pc2*qb2^2 - 2*pa2*pc3*qa1*qc2 - 2*pa2*pc3*qa2*qc1 + 2*pa2*pc3*qb1*qb2 + pa3^2*qc1^2 - pa3*pb1*qb1*qc3 - pa3*pb1*qb2*qc2 - pa3*pb1*qb3*qc1 - pa3*pb2*qb1*qc2 - pa3*pb2*qb2*qc1 - pa3*pb3*qb1*qc1 - 2*pa3*pc1*qa1*qc3 - 2*pa3*pc1*qa2*qc2 - 2*pa3*pc1*qa3*qc1 + 2*pa3*pc1*qb1*qb3 + pa3*pc1*qb2^2 - 2*pa3*pc2*qa1*qc2 - 2*pa3*pc2*qa2*qc1 + 2*pa3*pc2*qb1*qb2 - 2*pa3*pc3*qa1*qc1 + pa3*pc3*qb1^2 + pb1^2*qa3*qc3 + 2*pb1*pb2*qa2*qc3 + 2*pb1*pb2*qa3*qc2 + 2*pb1*pb3*qa1*qc3 + 2*pb1*pb3*qa2*qc2 + 2*pb1*pb3*qa3*qc1 - pb1*pc1*qa3*qb3 - pb1*pc2*qa2*qb3 - pb1*pc2*qa3*qb2 - pb1*pc3*qa1*qb3 - pb1*pc3*qa2*qb2 - pb1*pc3*qa3*qb1 + pb2^2*qa1*qc3 + pb2^2*qa2*qc2 + pb2^2*qa3*qc1 + 2*pb2*pb3*qa1*qc2 + 2*pb2*pb3*qa2*qc1 - pb2*pc1*qa2*qb3 - pb2*pc1*qa3*qb2 - pb2*pc2*qa1*qb3 - pb2*pc2*qa2*qb2 - pb2*pc2*qa3*qb1 - pb2*pc3*qa1*qb2 - pb2*pc3*qa2*qb1 + pb3^2*qa1*qc1 - pb3*pc1*qa1*qb3 - pb3*pc1*qa2*qb2 - pb3*pc1*qa3*qb1 - pb3*pc2*qa1*qb2 - pb3*pc2*qa2*qb1 - pb3*pc3*qa1*qb1 + pc1^2*qa3^2 + 4*pc1*pc2*qa2*qa3 + 4*pc1*pc3*qa1*qa3 + 2*pc1*pc3*qa2^2 + 2*pc2^2*qa1*qa3 + pc2^2*qa2^2 + 4*pc2*pc3*qa1*qa2 + pc3^2*qa1^2;
C(4) = 2*qc3*pa1^2*qc2 + 4*qc3*pa1*pa2*qc1 + 2*pa1*pa2*qc2^2 - qc3*pa1*pb1*qb2 - qb3*pa1*pb1*qc2 - qc3*pa1*pb2*qb1 - pa1*pb2*qb2*qc2 - qb3*pa1*pb2*qc1 - 2*qc3*pa1*pc1*qa2 + 2*qb3*pa1*pc1*qb2 - 2*qa3*pa1*pc1*qc2 - 2*qc3*pa1*pc2*qa1 - 2*pa1*pc2*qa2*qc2 + 2*qb3*pa1*pc2*qb1 + pa1*pc2*qb2^2 - 2*qa3*pa1*pc2*qc1 - 2*pc3*pa1*qa1*qc2 - 2*pc3*pa1*qa2*qc1 + 2*pc3*pa1*qb1*qb2 - pb3*pa1*qb1*qc2 - pb3*pa1*qb2*qc1 + 4*pa3*pa1*qc1*qc2 + 2*pa2^2*qc1*qc2 - qc3*pa2*pb1*qb1 - pa2*pb1*qb2*qc2 - qb3*pa2*pb1*qc1 - pa2*pb2*qb1*qc2 - pa2*pb2*qb2*qc1 - 2*qc3*pa2*pc1*qa1 - 2*pa2*pc1*qa2*qc2 + 2*qb3*pa2*pc1*qb1 + pa2*pc1*qb2^2 - 2*qa3*pa2*pc1*qc1 - 2*pa2*pc2*qa1*qc2 - 2*pa2*pc2*qa2*qc1 + 2*pa2*pc2*qb1*qb2 - 2*pc3*pa2*qa1*qc1 + pc3*pa2*qb1^2 - pb3*pa2*qb1*qc1 + 2*pa3*pa2*qc1^2 + qc3*pb1^2*qa2 + qa3*pb1^2*qc2 + 2*qc3*pb1*pb2*qa1 + 2*pb1*pb2*qa2*qc2 + 2*qa3*pb1*pb2*qc1 - qb3*pb1*pc1*qa2 - qa3*pb1*pc1*qb2 - qb3*pb1*pc2*qa1 - pb1*pc2*qa2*qb2 - qa3*pb1*pc2*qb1 - pc3*pb1*qa1*qb2 + 2*pb3*pb1*qa1*qc2 - pc3*pb1*qa2*qb1 + 2*pb3*pb1*qa2*qc1 - pa3*pb1*qb1*qc2 - pa3*pb1*qb2*qc1 + pb2^2*qa1*qc2 + pb2^2*qa2*qc1 - qb3*pb2*pc1*qa1 - pb2*pc1*qa2*qb2 - qa3*pb2*pc1*qb1 - pb2*pc2*qa1*qb2 - pb2*pc2*qa2*qb1 - pc3*pb2*qa1*qb1 + 2*pb3*pb2*qa1*qc1 - pa3*pb2*qb1*qc1 + 2*qa3*pc1^2*qa2 + 4*qa3*pc1*pc2*qa1 + 2*pc1*pc2*qa2^2 + 4*pc3*pc1*qa1*qa2 - pb3*pc1*qa1*qb2 - 2*pa3*pc1*qa1*qc2 - pb3*pc1*qa2*qb1 - 2*pa3*pc1*qa2*qc1 + 2*pa3*pc1*qb1*qb2 + 2*pc2^2*qa1*qa2 + 2*pc3*pc2*qa1^2 - pb3*pc2*qa1*qb1 - 2*pa3*pc2*qa1*qc1 + pa3*pc2*qb1^2;
C(3) = 2*qc3*pa1^2*qc1 + pa1^2*qc2^2 + 4*pa1*pa2*qc1*qc2 - qc3*pa1*pb1*qb1 - pa1*pb1*qb2*qc2 - qb3*pa1*pb1*qc1 - pa1*pb2*qb1*qc2 - pa1*pb2*qb2*qc1 - 2*qc3*pa1*pc1*qa1 - 2*pa1*pc1*qa2*qc2 + 2*qb3*pa1*pc1*qb1 + pa1*pc1*qb2^2 - 2*qa3*pa1*pc1*qc1 - 2*pa1*pc2*qa1*qc2 - 2*pa1*pc2*qa2*qc1 + 2*pa1*pc2*qb1*qb2 - 2*pc3*pa1*qa1*qc1 + pc3*pa1*qb1^2 - pb3*pa1*qb1*qc1 + 2*pa3*pa1*qc1^2 + pa2^2*qc1^2 - pa2*pb1*qb1*qc2 - pa2*pb1*qb2*qc1 - pa2*pb2*qb1*qc1 - 2*pa2*pc1*qa1*qc2 - 2*pa2*pc1*qa2*qc1 + 2*pa2*pc1*qb1*qb2 - 2*pa2*pc2*qa1*qc1 + pa2*pc2*qb1^2 + qc3*pb1^2*qa1 + pb1^2*qa2*qc2 + qa3*pb1^2*qc1 + 2*pb1*pb2*qa1*qc2 + 2*pb1*pb2*qa2*qc1 - qb3*pb1*pc1*qa1 - pb1*pc1*qa2*qb2 - qa3*pb1*pc1*qb1 - pb1*pc2*qa1*qb2 - pb1*pc2*qa2*qb1 - pc3*pb1*qa1*qb1 + 2*pb3*pb1*qa1*qc1 - pa3*pb1*qb1*qc1 + pb2^2*qa1*qc1 - pb2*pc1*qa1*qb2 - pb2*pc1*qa2*qb1 - pb2*pc2*qa1*qb1 + 2*qa3*pc1^2*qa1 + pc1^2*qa2^2 + 4*pc1*pc2*qa1*qa2 + 2*pc3*pc1*qa1^2 - pb3*pc1*qa1*qb1 - 2*pa3*pc1*qa1*qc1 + pa3*pc1*qb1^2 + pc2^2*qa1^2;
C(2) = 2*qc2*pa1^2*qc1 - qc2*pa1*pb1*qb1 - qb2*pa1*pb1*qc1 - 2*qc2*pa1*pc1*qa1 + 2*qb2*pa1*pc1*qb1 - 2*qa2*pa1*pc1*qc1 - 2*pc2*pa1*qa1*qc1 + pc2*pa1*qb1^2 - pb2*pa1*qb1*qc1 + 2*pa2*pa1*qc1^2 + qc2*pb1^2*qa1 + qa2*pb1^2*qc1 - qb2*pb1*pc1*qa1 - qa2*pb1*pc1*qb1 - pc2*pb1*qa1*qb1 + 2*pb2*pb1*qa1*qc1 - pa2*pb1*qb1*qc1 + 2*qa2*pc1^2*qa1 + 2*pc2*pc1*qa1^2 - pb2*pc1*qa1*qb1 - 2*pa2*pc1*qa1*qc1 + pa2*pc1*qb1^2;
C(1) = pa1^2*qc1^2 - pa1*pb1*qb1*qc1 - 2*pa1*pc1*qa1*qc1 + pa1*pc1*qb1^2 + pb1^2*qa1*qc1 - pb1*pc1*qa1*qb1 + pc1^2*qa1^2;

x1 = roots(C);

% x2 can be found by solving a linear equation
% careful, this may not be numerically stable
x2 = NaN(size(x1));
N = length(x2);

for i = 1:N
    Pa_i = Pa(1)*x1(i)^2 + Pa(2)*x1(i) + Pa(3);
    Pb_i = Pb(1)*x1(i)^2 + Pb(2)*x1(i) + Pb(3);
    Pc_i = Pc(1)*x1(i)^2 + Pc(2)*x1(i) + Pc(3);
    Qa_i = Qa(1)*x1(i)^2 + Qa(2)*x1(i) + Qa(3);
    Qb_i = Qb(1)*x1(i)^2 + Qb(2)*x1(i) + Qb(3);
    Qc_i = Qc(1)*x1(i)^2 + Qc(2)*x1(i) + Qc(3);

    x2(i) = -(Pc_i/Pa_i - Qc_i/Qa_i) / (Pb_i/Pa_i - Qb_i/Qa_i);
end

end