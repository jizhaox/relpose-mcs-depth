function generate_construct_coef_matrix()
%% generate the coefficients of the equations using symbolic computation

% unknown variables
syms lmd1 lmd2 s real
syms x y z real

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

x11 = [x11_1; x11_2; x11_3];
x12 = [x12_1; x12_2; x12_3];
x21 = [x21_1; x21_2; x21_3];
x22 = [x22_1; x22_2; x22_3];
A1 = [a1_11, a1_12, 0; a1_21, a1_22, 0; 0, 0, 0];
A2 = [a2_11, a2_12, 0; a2_21, a2_22, 0; 0, 0, 0];
Q11 = [q11_11, q11_12, q11_13; q11_21, q11_22, q11_23; q11_31, q11_32, q11_33];
Q12 = [q12_11, q12_12, q12_13; q12_21, q12_22, q12_23; q12_31, q12_32, q12_33];
s11 = [s11_1; s11_2; s11_3];
s12 = [s12_1; s12_2; s12_3]*s; % unknown scale
Q21 = [q21_11, q21_12, q21_13; q21_21, q21_22, q21_23; q21_31, q21_32, q21_33];
Q22 = [q22_11, q22_12, q22_13; q22_21, q22_22, q22_23; q22_31, q22_32, q22_33];
s21 = [s21_1; s21_2; s21_3];
s22 = [s22_1; s22_2; s22_3]*s; % unknown scale
R = [1+x^2-y^2-z^2, 2*x*y-2*z, 2*y+2*x*z;
     2*x*y+2*z, 1-x^2+y^2-z^2, 2*y*z-2*x;
     2*x*z-2*y, 2*x+2*y*z, 1-x^2-y^2+z^2];

p11 = [p11_1; p11_2; p11_3];
p12 = [p12_1; p12_2; p12_3];
q11 = cross(s11, p11);
q12 = cross(s12, p12);
pq11 = cross(p11, q11);
pq12 = cross(p12, q12);

t1 = pq11 + lmd1 * p11;
t2 = pq12 + lmd2 * p12;

E1 = Q12'*(R*skew(s11-t1) + skew(t2-s12)*R)*Q11;
%eq_epi_1 = x12'*E1*x11;
eq_aff_1 = E1'*x12 + A1'*E1*x11;
eq_aff_1 = eq_aff_1(1:2);

E2 = Q22'*(R*skew(s21-t1) + skew(t2-s22)*R)*Q21;
eq_epi_2 = x22'*E2*x21;
eq_aff_2 = E2'*x22 + A2'*E2*x21;
eq_aff_2 = eq_aff_2(1:2);

eqs = [eq_aff_1;eq_epi_2;eq_aff_2];
[A_sym, b] = equationsToMatrix(eqs, [lmd1, lmd2, s]);
A = [A_sym, -b];
% check
simplify(A*[lmd1; lmd2; s; 1] - eqs)

filename = 'str1.m';
fid = fopen(filename, 'r');
text1 = fread(fid, '*char')';
text1 = strrep(text1, sprintf('\r'), '');
fclose(fid);
fid = fopen('construct_coef_epi.m', 'wt');
fprintf(fid, '%s\n', 'function C = construct_coef_epi(p11, p12, s11, s12, x21, x22, Q21, s21, Q22, s22)');
fprintf(fid, '%s\n\n', text1);
fprintf(fid, '%s\n', 'C = zeros(1, 4, 10);');
cnt = 1;
for i = 3:3
    for j = 1:size(A,2)
        [coef_tmp, term_tmp] = coeffs(A(i,j), [x,y,z]);
%         disp(term_tmp)
         for k = 1:10
             fprintf(fid, '%s\n', ['C(' num2str(cnt) ', ' num2str(j) ', ' num2str(k) ') = ' char(coef_tmp(k)) ';']);
         end
    end
    cnt = cnt + 1;
end
fclose(fid);

filename = 'str2.m';
fid = fopen(filename, 'r');
text2 = fread(fid, '*char')';
text2 = strrep(text2, sprintf('\r'), '');
fclose(fid);
fid = fopen('construct_coef_aff.m', 'wt');
fprintf(fid, '%s\n', 'function C = construct_coef_aff(p11, p12, s11, s12, x21, x22, A2, Q21, s21, Q22, s22)');
fprintf(fid, '%s\n\n', text2);
fprintf(fid, '%s\n', 'C = zeros(2, 4, 10);');
cnt = 1;
for i = 4:5
    for j = 1:size(A,2)
        [coef_tmp, term_tmp] = coeffs(A(i,j), [x,y,z]);
%         disp(term_tmp)
         for k = 1:10
             fprintf(fid, '%s\n', ['C(' num2str(cnt) ', ' num2str(j) ', ' num2str(k) ') = ' char(coef_tmp(k)) ';']);
         end
    end
    cnt = cnt + 1;
end
fclose(fid);