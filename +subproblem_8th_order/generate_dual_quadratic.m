syms pa pb pc qa qb qc [1 3]
syms x1 x2

Pa = pa1*x1^2 + pa2*x1 + pa3;
Pb = pb1*x1^2 + pb2*x1 + pb3;
Pc = pc1*x1^2 + pc2*x1 + pc3;
Qa = qa1*x1^2 + qa2*x1 + qa3;
Qb = qb1*x1^2 + qb2*x1 + qb3;
Qc = qc1*x1^2 + qc2*x1 + qc3;

sylvester = [Pa Pb Pc 0
            0  Pa Pb Pc
            Qa Qb Qc 0
            0  Qa Qb Qc];

S = det(sylvester);
C = coeffs(S, x1);