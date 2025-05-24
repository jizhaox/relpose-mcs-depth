function generate_construct_coeff()
syms x1_1 x1_2 x1_3 real
syms x2_1 x2_2 x2_3 real
syms a11 a12 a21 a22 real
syms q1_11 q1_12 q1_13 q1_21 q1_22 q1_23 q1_31 q1_32 q1_33 real
syms q2_11 q2_12 q2_13 q2_21 q2_22 q2_23 q2_31 q2_32 q2_33 real
syms s1_1 s1_2 s1_3 real
syms s2_1 s2_2 s2_3 real
syms lmd1 lmd2 real
syms qy real
syms p1_1 p1_2 p1_3 real
syms p2_1 p2_2 p2_3 real
syms qp1_1 qp1_2 qp1_3 real
syms qp2_1 qp2_2 qp2_3 real

x1 = [x1_1; x1_2; x1_3];
x2 = [x2_1; x2_2; x2_3];
A = [a11, a12, 0; a21, a22, 0; 0, 0, 0];
Q1 = [q1_11, q1_12, q1_13; q1_21, q1_22, q1_23; q1_31, q1_32, q1_33];
Q2 = [q2_11, q2_12, q2_13; q2_21, q2_22, q2_23; q2_31, q2_32, q2_33];
s1 = [s1_1; s1_2; s1_3];
s2 = [s2_1; s2_2; s2_3];
R = [1-qy^2, 0, -2*qy; 0, 1+qy^2, 0; 2*qy, 0, 1-qy^2];
p1 = [p1_1; p1_2; p1_3];
p2 = [p2_1; p2_2; p2_3];
qp1 = [qp1_1; qp1_2; qp1_3];
qp2 = [qp2_1; qp2_2; qp2_3];
t1 = qp1 + lmd1 * p1;
t2 = qp2 + lmd2 * p2;

E2 = Q2'*(R*skew(s1-t1) + skew(t2-s2)*R)*Q1;
eq_epi = x2'*E2*x1;
eq_aff = E2'*x2 + A'*E2*x1;
eqs = [eq_epi; eq_aff(1:2)];

[A_sym, b] = equationsToMatrix(eqs, [lmd1, lmd2]);
A = [A_sym, -b];
% check
%simplify(A*[lmd1; lmd2; 1] - eqs)

syms c11_2 c11_1 c11_0 real
syms c12_2 c12_1 c12_0 real
syms c13_2 c13_1 c13_0 real
syms c21_2 c21_1 c21_0 real
syms c22_2 c22_1 c22_0 real
syms c23_2 c23_1 c23_0 real
syms c31_2 c31_1 c31_0 real
syms c32_2 c32_1 c32_0 real
syms c33_2 c33_1 c33_0 real

c11 = c11_2*qy^2 + c11_1*qy + c11_0;
c12 = c12_2*qy^2 + c12_1*qy + c12_0;
c13 = c13_2*qy^2 + c13_1*qy + c13_0;
c21 = c21_2*qy^2 + c21_1*qy + c21_0;
c22 = c22_2*qy^2 + c22_1*qy + c22_0;
c23 = c23_2*qy^2 + c23_1*qy + c23_0;
c31 = c31_2*qy^2 + c31_1*qy + c31_0;
c32 = c32_2*qy^2 + c32_1*qy + c32_0;
c33 = c33_2*qy^2 + c33_1*qy + c33_0;

F = [c11, c12, c13;
    c21, c22, c23;
    c31, c32, c33];
[coef_eq, term_eq] = coeffs(det(F),qy);

% save
fid = fopen('construct_coeff.m', 'wt');
for i = 1:3
    for j = 1:3
        [coef_tmp, term_tmp] = coeffs(A(i,j), qy);
        for k = 1:3
            fprintf(fid, '%s\n', ['c' num2str(i) num2str(j) '_' num2str(3-k) ' = ' char(coef_tmp(k)) ';']);
        end
    end
end
fprintf(fid, '%s\n', ['coef_eq = ' char(coef_eq) ';']);
fclose(fid);