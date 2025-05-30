function generate_construct_coeff()
%% generate the coefficients of the equations using symbolic computation

% xij: image coordinates of i-th AC at j-th view
syms x11_1 x11_2 x11_3 real
syms x12_1 x12_2 x12_3 real
syms x21_1 x21_2 x21_3 real
syms x22_1 x22_2 x22_3 real
% Ai: affine matrix of i-th AC
syms a1_11 a1_12 a1_21 a1_22 real
syms a2_11 a2_12 a2_21 a2_22 real
% Qij, sij: extrinsic parameters of i-th AC at j-th view
syms q11_11 q11_12 q11_13 q11_21 q11_22 q11_23 q11_31 q11_32 q11_33 real
syms q12_11 q12_12 q12_13 q12_21 q12_22 q12_23 q12_31 q12_32 q12_33 real
syms s11_1 s11_2 s11_3 real
syms s12_1 s12_2 s12_3 real
syms q21_11 q21_12 q21_13 q21_21 q21_22 q21_23 q21_31 q21_32 q21_33 real
syms q22_11 q22_12 q22_13 q22_21 q22_22 q22_23 q22_31 q22_32 q22_33 real
syms s21_1 s21_2 s21_3 real
syms s22_1 s22_2 s22_3 real
% pij, qij: plucker line of i-th AC at j-th view
syms p11_1 p11_2 p11_3 real
syms p12_1 p12_2 p12_3 real
syms pq11_1 pq11_2 pq11_3 real
syms pq12_1 pq12_2 pq12_3 real

% unknown variables
syms lmd1 lmd2 real
syms qy real

x11 = [x11_1; x11_2; x11_3];
x12 = [x12_1; x12_2; x12_3];
x21 = [x21_1; x21_2; x21_3];
x22 = [x22_1; x22_2; x22_3];
A1 = [a1_11, a1_12, 0; a1_21, a1_22, 0; 0, 0, 0];
A2 = [a2_11, a2_12, 0; a2_21, a2_22, 0; 0, 0, 0];
Q11 = [q11_11, q11_12, q11_13; q11_21, q11_22, q11_23; q11_31, q11_32, q11_33];
Q12 = [q12_11, q12_12, q12_13; q12_21, q12_22, q12_23; q12_31, q12_32, q12_33];
s11 = [s11_1; s11_2; s11_3];
s12 = [s12_1; s12_2; s12_3];
Q21 = [q21_11, q21_12, q21_13; q21_21, q21_22, q21_23; q21_31, q21_32, q21_33];
Q22 = [q22_11, q22_12, q22_13; q22_21, q22_22, q22_23; q22_31, q22_32, q22_33];
s21 = [s21_1; s21_2; s21_3];
s22 = [s22_1; s22_2; s22_3];
R = [1-qy^2, 0, -2*qy; 0, 1+qy^2, 0; 2*qy, 0, 1-qy^2];
p11 = [p11_1; p11_2; p11_3];
p12 = [p12_1; p12_2; p12_3];
pq11 = [pq11_1; pq11_2; pq11_3];
pq12 = [pq12_1; pq12_2; pq12_3];
t1 = pq11 + lmd1 * p11;
t2 = pq12 + lmd2 * p12;

E1 = Q12'*(R*skew(s11-t1) + skew(t2-s12)*R)*Q11;
%eq_epi_1 = x12'*E1*x11;
eq_aff_1 = E1'*x12 + A1'*E1*x11;

E2 = Q22'*(R*skew(s21-t1) + skew(t2-s22)*R)*Q21;
eq_epi_2 = x22'*E2*x21;
eq_aff_2 = E2'*x22 + A2'*E2*x21;

eqs = [eq_aff_1(1); eq_epi_2; eq_aff_2(1)];

[A_sym, b] = equationsToMatrix(eqs, [lmd1, lmd2]);
A = [A_sym, -b];
% check
simplify(A*[lmd1; lmd2; 1] - eqs)

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
            fprintf(fid, '%s\n', ['c' num2str(i) num2str(j) '_' num2str(3-k) ' = ' char(simplify(coef_tmp(k))) ';']);
        end
    end
end
fprintf(fid, '%s\n', ['coef_eq = ' char(simplify(coef_eq)) ';']);
fclose(fid);