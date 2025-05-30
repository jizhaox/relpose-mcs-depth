function [cay_sols, t_sols, R_sols, s_sols, cay_sols_all] = solver_depth_3ac(Image1, Image2, At, R_cam_first, R_cam_second, t_cam_first, t_cam_second)

cay_sols = [];
t_sols = [];
R_sols = [];
s_sols = [];

% xij: image coordinates of i-th AC at j-th view
x11 = Image1(:, 1);
x12 = Image2(:, 1);
x21 = Image1(:, 2);
x22 = Image2(:, 2);
x31 = Image1(:, 3);
x32 = Image2(:, 3);

% Ai: affine matrix of i-th AC
A1 = At(:, :, 1);
A2 = At(:, :, 2);
% A3 = At(:, :, 3); % redundancy

% Qij, sij: extrinsic parameters of i-th AC at j-th view
Q11 = R_cam_first(:, :, 1);
s11 = t_cam_first(:, 1);
Q12 = R_cam_second(:, :, 1);
s12 = t_cam_second(:, 1);

Q21 = R_cam_first(:, :, 2);
s21 = t_cam_first(:, 2);
Q22 = R_cam_second(:, :, 2);
s22 = t_cam_second(:, 2);

Q31 = R_cam_first(:, :, 3);
s31 = t_cam_first(:, 3);
Q32 = R_cam_second(:, :, 3);
s32 = t_cam_second(:, 3);

% pij, qij: plucker line of i-th AC at j-th view
% pij, qij: plucker line of i-th AC at j-th view
p11 = Q11*x11;
p11 = p11/norm(p11(:));
p12 = Q12*x12;
p12 = p12/norm(p12(:));

p21 = Q21*x21;
p21 = p21/norm(p21(:));
p22 = Q22*x22;
p22 = p22/norm(p22(:));

p31 = Q31*x31;
p31 = p31/norm(p31(:));
p32 = Q32*x32;
p32 = p32/norm(p32(:));

%% choose the 1st AC to build the reference
% The fist 4 inputs are about the selected AC, the other inputs are about 
% the AC that used to form the coefficient matrix
M1 = construct_coef_aff(p11, p12, s11, s12, ...
    x11, x12, A1, Q11, s11, Q12, s12);
M2 = construct_coef_epi(p11, p12, s11, s12, ...
    x21, x22, Q21, s21, Q22, s22);
M3 = construct_coef_aff(p11, p12, s11, s12, ...
    x21, x22, A2, Q21, s21, Q22, s22);
M4 = construct_coef_epi(p11, p12, s11, s12, ...
    x31, x32, Q31, s31, Q32, s32);
F1 = cat(1, M1, M2, M3, M4);
%% choose the 2nd AC to build the reference
% The fist 4 inputs are about the selected AC, the other inputs are about 
% the AC that used to form the coefficient matrix
M1 = construct_coef_aff(p21, p22, s21, s22, ...
    x21, x22, A2, Q21, s21, Q22, s22);
M2 = construct_coef_epi(p21, p22, s21, s22, ...
    x11, x12, Q11, s11, Q12, s12);
M3 = construct_coef_aff(p21, p22, s21, s22, ...
    x11, x12, A1, Q11, s11, Q12, s12);
M4 = construct_coef_epi(p21, p22, s21, s22, ...
    x31, x32, Q31, s31, Q32, s32);
F2 = cat(1, M1, M2, M3, M4);
%% choose the 3rd AC to build the reference
% The fist 4 inputs are about the selected AC, the other inputs are about 
% the AC that used to form the coefficient matrix
M1 = construct_coef_epi(p31, p32, s31, s32, ...
    x11, x12, Q11, s11, Q12, s12);
M2 = construct_coef_aff(p31, p32, s31, s32, ...
    x11, x12, A1, Q11, s11, Q12, s12);
M3 = construct_coef_epi(p31, p32, s31, s32, ...
    x21, x22, Q21, s21, Q22, s22);
M4 = construct_coef_aff(p31, p32, s31, s32, ...
    x21, x22, A2, Q21, s21, Q22, s22);
F3 = cat(1, M1, M2, M3, M4);

%%
idx1_4by4 = [
     1     2     3     4     1     2     3     4
     1     2     3     5     1     2     3     4
     1     2     4     5     1     2     3     4
     1     3     4     5     1     2     3     4
     1     3     4     6     1     2     3     4
     1     3     5     6     1     2     3     4
     1     4     5     6     1     2     3     4
     2     3     4     5     1     2     3     4
     2     3     4     6     1     2     3     4
     2     3     5     6     1     2     3     4
     2     4     5     6     1     2     3     4
     3     4     5     6     1     2     3     4];
idx1_2by2 = [
     1     2     1     2
     1     2     1     3
     1     2     2     4
     1     2     3     4
     1     6     1     2
     1     6     1     3
     1     6     2     4
     1     6     3     4
     2     6     1     2
     2     6     1     3
     2     6     2     4
     2     6     3     4];
idx2_4by4 = [
     1     2     3     4     1     2     3     4
     1     2     3     5     1     2     3     4
     1     2     3     6     1     2     3     4
     1     2     4     5     1     2     3     4
     1     2     4     6     1     2     3     4
     1     2     5     6     1     2     3     4
     1     3     4     5     1     2     3     4
     1     3     4     6     1     2     3     4
     1     3     5     6     1     2     3     4
     1     4     5     6     1     2     3     4
     2     3     4     5     1     2     3     4
     2     3     4     6     1     2     3     4
     2     3     5     6     1     2     3     4
     2     4     5     6     1     2     3     4];
idx2_2by2 = [
     1     2     1     2
     1     2     1     3
     1     2     2     4
     1     2     3     4];
idx3_4by4 = [
     1     2     4     5     1     2     3     4
     1     2     4     6     1     2     3     4
     1     2     5     6     1     2     3     4
     1     3     4     5     1     2     3     4
     1     3     4     6     1     2     3     4
     1     3     5     6     1     2     3     4
     1     4     5     6     1     2     3     4
     2     3     4     5     1     2     3     4
     2     3     4     6     1     2     3     4
     2     3     5     6     1     2     3     4
     2     4     5     6     1     2     3     4
     3     4     5     6     1     2     3     4];
idx3_2by2 = [
     1     2     1     2
     1     2     1     3
     1     2     2     4
     1     2     3     4
     1     3     1     2
     1     3     1     3
     1     3     2     4
     1     3     3     4
     2     3     1     2
     2     3     1     3
     2     3     2     4
     2     3     3     4];

% equations from 4*4 determinants
cnt = 1;
C1 = zeros(165, size(idx1_4by4, 1)+size(idx2_4by4, 1)+size(idx3_4by4, 1));
for i = 1:size(idx1_4by4, 1)
    %
    ir = idx1_4by4(i, 1:4); % index for row
    ic = idx1_4by4(i, 5:8); % index for column
    C1(:, cnt) = det4x4poly(F1(ir, ic, :));
    cnt = cnt + 1;
end
for i = 1:size(idx2_4by4, 1)
    %
    ir = idx2_4by4(i, 1:4); % index for row
    ic = idx2_4by4(i, 5:8); % index for column
    C1(:, cnt) = det4x4poly(F2(ir, ic, :));
    cnt = cnt + 1;
end
for i = 1:size(idx3_4by4, 1)
    %
    ir = idx3_4by4(i, 1:4); % index for row
    ic = idx3_4by4(i, 5:8); % index for column
    C1(:, cnt) = det4x4poly(F3(ir, ic, :));
    cnt = cnt + 1;
end
C1 = C1(1:end-1, :);

% equations from 2*2 determinants
cnt = 1;
C2 = zeros(35, size(idx1_2by2, 1)+size(idx2_2by2, 1)+size(idx3_2by2, 1));
for i = 1:size(idx1_2by2, 1)
    %
    ir = idx1_2by2(i, 1:2); % index for row
    ic = idx1_2by2(i, 3:4); % index for column
    C2(:, cnt) = det2x2poly(F1(ir, ic, :));
    cnt = cnt + 1;
end
for i = 1:size(idx2_2by2, 1)
    %
    ir = idx2_2by2(i, 1:2); % index for row
    ic = idx2_2by2(i, 3:4); % index for column
    C2(:, cnt) = det2x2poly(F2(ir, ic, :));
    cnt = cnt + 1;
end
for i = 1:size(idx3_2by2, 1)
    %
    ir = idx3_2by2(i, 1:2); % index for row
    ic = idx3_2by2(i, 3:4); % index for column
    C2(:, cnt) = det2x2poly(F3(ir, ic, :));
    cnt = cnt + 1;
end

cay_sols_all = solver_depth_3ac_core([C1(:); C2(:)]);
idx = find(all(imag(cay_sols_all) == 0, 1));
if isempty(idx)
    return;
end

n = numel(idx);
cay_sols = cay_sols_all(:, idx);
t_sols = zeros(3, n);
R_sols = zeros(3, 3, n);
s_sols = zeros(n, 1);
% estimate scale and translation
for k = 1:n
    cay = cay_sols(:, k);
    x = cay(1);
    y = cay(2);
    z = cay(3);
    mom = [x^2, x*y, x*z, x, y^2, y*z, y, z^2, z, 1];
    F = zeros(size(F1,1), size(F1,2));
    for i = 1:size(F1,1)
        for j = 1:size(F1,2)
            F(i,j) = mom*squeeze(F1(i,j,:));
        end
    end
    [U, S, V] = svd(F, 'econ');
    [~, idx] = min(diag(S));
    x = V(:, idx);
    x = x / x(end);
    s_sols(k) = x(3);

    lmd = x(1:2);
    [p11, ~, pq11] = plucker_line(x11, Q11, s11);
    [p12, ~, pq12] = plucker_line(x12, Q12, s12*s_sols(k));
    t1 = pq11 + lmd(1) * p11;
    t2 = pq12 + lmd(2) * p12;

    R = cayley_rotation(cay);
    R_sols(:, :, k) = R;
    t_sols(:, k) = -R*t1 + t2;
end


%% check
% x = cay_gt(1);
% y = cay_gt(2);
% z = cay_gt(3);
% mm = [x^2, x*y, x*z, x, y^2, y*z, y, z^2, z, 1];
% monom_str = {'x^8','x^7*y','x^6*y^2','x^5*y^3','x^4*y^4','x^3*y^5','x^2*y^6','x*y^7','y^8','x^7*z','x^6*y*z','x^5*y^2*z','x^4*y^3*z','x^3*y^4*z','x^2*y^5*z','x*y^6*z','y^7*z','x^6*z^2','x^5*y*z^2','x^4*y^2*z^2','x^3*y^3*z^2','x^2*y^4*z^2','x*y^5*z^2','y^6*z^2','x^5*z^3','x^4*y*z^3','x^3*y^2*z^3','x^2*y^3*z^3','x*y^4*z^3','y^5*z^3','x^4*z^4','x^3*y*z^4','x^2*y^2*z^4','x*y^3*z^4','y^4*z^4','x^3*z^5','x^2*y*z^5','x*y^2*z^5','y^3*z^5','x^2*z^6','x*y*z^6','y^2*z^6','x*z^7','y*z^7','z^8','x^7','x^6*y','x^5*y^2','x^4*y^3','x^3*y^4','x^2*y^5','x*y^6','y^7','x^6*z','x^5*y*z','x^4*y^2*z','x^3*y^3*z','x^2*y^4*z','x*y^5*z','y^6*z','x^5*z^2','x^4*y*z^2','x^3*y^2*z^2','x^2*y^3*z^2','x*y^4*z^2','y^5*z^2','x^4*z^3','x^3*y*z^3','x^2*y^2*z^3','x*y^3*z^3','y^4*z^3','x^3*z^4','x^2*y*z^4','x*y^2*z^4','y^3*z^4','x^2*z^5','x*y*z^5','y^2*z^5','x*z^6','y*z^6','z^7','x^6','x^5*y','x^4*y^2','x^3*y^3','x^2*y^4','x*y^5','y^6','x^5*z','x^4*y*z','x^3*y^2*z','x^2*y^3*z','x*y^4*z','y^5*z','x^4*z^2','x^3*y*z^2','x^2*y^2*z^2','x*y^3*z^2','y^4*z^2','x^3*z^3','x^2*y*z^3','x*y^2*z^3','y^3*z^3','x^2*z^4','x*y*z^4','y^2*z^4','x*z^5','y*z^5','z^6','x^5','x^4*y','x^3*y^2','x^2*y^3','x*y^4','y^5','x^4*z','x^3*y*z','x^2*y^2*z','x*y^3*z','y^4*z','x^3*z^2','x^2*y*z^2','x*y^2*z^2','y^3*z^2','x^2*z^3','x*y*z^3','y^2*z^3','x*z^4','y*z^4','z^5','x^4','x^3*y','x^2*y^2','x*y^3','y^4','x^3*z','x^2*y*z','x*y^2*z','y^3*z','x^2*z^2','x*y*z^2','y^2*z^2','x*z^3','y*z^3','z^4','x^3','x^2*y','x*y^2','y^3','x^2*z','x*y*z','y^2*z','x*z^2','y*z^2','z^3','x^2','x*y','y^2','x*z','y*z','z^2','x','y','z','1'};
% monoms = str2sym(monom_str);
% vals = subs(monoms, {'x', 'y', 'z'}, {x, y, z});
% vals_num = double(vals); 
% vals_num(1:end-1)*C1
% vals_num(131:end)*C2

