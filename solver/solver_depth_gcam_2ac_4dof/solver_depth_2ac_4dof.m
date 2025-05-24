function [R_sol, t_sol, theta_y_sol] = solver_depth_2ac_4dof(Image1, Image2, At, R_cam, t_cam, match_type)

%% generate the coefficients by symbolic computation
% notice: run only once
if 0
    generate_construct_coeff();
end

warning off

match_info = match_type_gcam_ac(match_type);
% choose arbitray one AC to build the reference
for ii = 1:1
    x1 = Image1(:, ii);
    x2 = Image2(:, ii);
    A = At(:, :, ii);
    % extract extrinsic parameters
    idx1 = match_info{ii}.idx1;
    idx2 = match_info{ii}.idx2;
    Q1 = R_cam(:, :, idx1);
    s1 = t_cam(:, idx1);
    Q2 = R_cam(:, :, idx2);
    s2 = t_cam(:, idx2);
    % line
    [p1, q1, qp1] = plucker_line(x1, Q1, s1);
    [p2, q2, qp2] = plucker_line(x2, Q2, s2);

    x1_1 = x1(1);
    x1_2 = x1(2);
    x1_3 = x1(3);
    x2_1 = x2(1);
    x2_2 = x2(2);
    x2_3 = x2(3);
    a11 = A(1,1);
    a12 = A(1,2);
    a21 = A(2,1);
    a22 = A(2,2);
    q1_11 = Q1(1,1);
    q1_12 = Q1(1,2);
    q1_13 = Q1(1,3);
    q1_21 = Q1(2,1);
    q1_22 = Q1(2,2);
    q1_23 = Q1(2,3);
    q1_31 = Q1(3,1);
    q1_32 = Q1(3,2);
    q1_33 = Q1(3,3);
    q2_11 = Q2(1,1);
    q2_12 = Q2(1,2);
    q2_13 = Q2(1,3);
    q2_21 = Q2(2,1);
    q2_22 = Q2(2,2);
    q2_23 = Q2(2,3);
    q2_31 = Q2(3,1);
    q2_32 = Q2(3,2);
    q2_33 = Q2(3,3);
    s1_1 = s1(1);
    s1_2 = s1(2);
    s1_3 = s1(3);
    s2_1 = s2(1);
    s2_2 = s2(2);
    s2_3 = s2(3);
    p1_1 = p1(1);
    p1_2 = p1(2);
    p1_3 = p1(3);
    p2_1 = p2(1);
    p2_2 = p2(2);
    p2_3 = p2(3);
    qp1_1 = qp1(1);
    qp1_2 = qp1(2);
    qp1_3 = qp1(3);
    qp2_1 = qp2(1);
    qp2_2 = qp2(2);
    qp2_3 = qp2(3);
    
    %
    construct_coeff;
    
    % solver qy
    qy_all = roots(coef_eq);
    idx = not(imag( qy_all ));
    qy_sol = qy_all(idx);
    theta_y_sol = atan(qy_sol)*2;

    n_sol = numel(qy_sol);
    R_sol = zeros(3, 3, n_sol);
    t_sol = zeros(3, n_sol);
    for k = 1:n_sol
        qy = qy_sol(k);
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
        lmd = -F(:, 1:2)\F(:, 3);
        t1 = qp1 + lmd(1) * p1;
        t2 = qp2 + lmd(2) * p2;
        
        R = [1-qy^2, 0, -2*qy; 0, 1+qy^2, 0; 2*qy, 0, 1-qy^2] / (1+qy^2);
        R_sol(:, :, k) = R;
        t_sol(:, k) = -R*t1 + t2;
    end

end

function [p, q, q_cross_p] = plucker_line(x, Q, s)
p = Q*x;
p = p/norm(p(:));
q = cross(s, p);
q_cross_p = cross(q, p);