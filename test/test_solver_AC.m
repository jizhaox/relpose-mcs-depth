clear;
addpath('../solver')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('--------------2AC method & inter-camera case--------------')
%% generate synthetic data
[Image1, Image2, At, R_cam, t_cam, R_gt, cay_gt, t_gt] = generate_2AC_data('inter');

%% check data
% For noise-free data, the residual of epipolar constraint and affine transformation constraints should be small (Eqs.(7) and (8) in the paper)
err_epipolar = zeros(2, 1);
err_affinetransformation = zeros(2, 2);
Tf_gt = [R_gt t_gt; 0 0 0 1];
Tf_cam{1} = [R_cam(:,:,1), t_cam(:,1); 0, 0, 0, 1];
Tf_cam{2} = [R_cam(:,:,2), t_cam(:,2); 0, 0, 0, 1];
for ii = 1:2
    if ii <= 1
        idx1 = 1;
        idx2 = 2;
    else
        idx1 = 2;
        idx2 = 1;
    end
    x1 = Image1(:, ii);
    x2 = Image2(:, ii);
    % relative pose between time i and time j 
    Hij = Tf_cam{idx2}\Tf_gt*(Tf_cam{idx1});
    Rij = Hij(1:3,1:3);
    tij = Hij(1:3,4);
    % epipolar constraint of PC
    x = tij;
    skew_tij=[0 -x(3) x(2); x(3) 0 -x(1); -x(2) x(1) 0];
    E =skew_tij*Rij;
    err_epipolar(ii) = x2'*E*x1;
    % affine transformation constraints
    At_33(1:2,1:2) = At(:,:,ii);
    At_33(3,3) = 1;
    equationerror = E'*x2 + At_33'*E*x1;
    err_affinetransformation(:,ii) = equationerror(1:2,1);
end
disp(['The maximum residual of epipolar constraint: ' num2str(max(abs(err_epipolar(:))))])
disp(['The maximum residual of affine transformation constraints: ' num2str(max(abs(err_affinetransformation(:))))])
    
%% run solver
% set the 6th parameter as 0 (default), 48 solutions
[cay_sols, t_sols, R_sols, cay_sols_all] = solver_depth_inter_2ac(Image1(1:2,:), Image2(1:2,:), At, R_cam, t_cam, 0);
cay_sol = find_solution(cay_sols, cay_gt);
t_sol = find_solution(t_sols, t_gt);
cay_sol, t_sol, cay_gt, t_gt

% set the 6th parameter as 1, 56 solutions.
% according to our evaluation, 56-solution configuration has better numerical stability
[cay_sols, t_sols, R_sols, cay_sols_all] = solver_depth_inter_2ac(Image1(1:2,:), Image2(1:2,:), At, R_cam, t_cam, 1);
cay_sol = find_solution(cay_sols, cay_gt);
t_sol = find_solution(t_sols, t_gt);
cay_sol, t_sol, cay_gt, t_gt


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('--------------2AC method & intra-camera case--------------')
[Image1, Image2, At, R_cam, t_cam, R_gt, cay_gt, t_gt] = generate_2AC_data('intra');

%% check data
% For noise-free data, the residual of epipolar constraint and affine transformation constraints should be small (Eqs.(7) and (8) in the paper)
err_epipolar = zeros(2, 1);
err_affinetransformation = zeros(2, 2);
Tf_gt = [R_gt t_gt; 0 0 0 1];
Tf_cam{1} = [R_cam(:,:,1), t_cam(:,1); 0, 0, 0, 1];
Tf_cam{2} = [R_cam(:,:,2), t_cam(:,2); 0, 0, 0, 1];
for ii = 1:2
    if ii <= 1
        idx1 = 1;
        idx2 = 1;
    else
        idx1 = 2;
        idx2 = 2;
    end
    x1 = Image1(:, ii);
    x2 = Image2(:, ii);
    % relative pose between time i and time j 
    Hij = Tf_cam{idx2}\Tf_gt*(Tf_cam{idx1});
    Rij = Hij(1:3,1:3);
    tij = Hij(1:3,4);
    % epipolar constraint of PC
    x = tij;
    skew_tij=[0 -x(3) x(2) ; x(3) 0 -x(1) ; -x(2) x(1) 0 ];
    E =skew_tij*Rij;
    err_epipolar(ii) = x2'*E*x1;
    % affine transformation constraints
    At_33(1:2,1:2) = At(:,:,ii);
    At_33(3,3) = 1;
    equationerror = E'*x2 + At_33'*E*x1;
    err_affinetransformation(:,ii) = equationerror(1:2,1); 
end
disp(['The maximum residual of epipolar geometry: ' num2str(max(abs(err_epipolar(:))))])
disp(['The maximum residual of affine transformation constraints: ' num2str(max(abs(err_affinetransformation(:))))])

%% run solver
% 48 solutions
[cay_sols, t_sols, R_sols, cay_sols_all] = solver_depth_intra_2ac(Image1(1:2,:), Image2(1:2,:), At, R_cam, t_cam);
cay_sol = find_solution(cay_sols, cay_gt);
t_sol = find_solution(t_sols, t_gt);
cay_sol, t_sol, cay_gt, t_gt

% add coordinate transformation
% according to our evaluation, it has better numerical stability for some cases
[cay_sols, t_sols, R_sols] = solver_depth_intra_2ac_enhance(Image1(1:2,:), Image2(1:2,:), At, R_cam, t_cam);
cay_sol = find_solution(cay_sols, cay_gt);
t_sol = find_solution(t_sols, t_gt);
cay_sol, t_sol, cay_gt, t_gt


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('--------------2AC method & single camera--------------')
n_point = 2;
[Image1, Image2, At, R_gt, cay_gt, t_gt] = generate_2AC_single_data(n_point);

%% check data
% For noise-free data, the residual of epipolar constraint and affine transformation constraints should be small (Eqs.(7) and (8) in the paper)
err_epipolar = zeros(2, 1);
err_affinetransformation = zeros(2, 2);

for ii = 1:2
    x1 = Image1(:, ii);
    x2 = Image2(:, ii);
    % epipolar constraint of PC
    x = t_gt;
    skew_tij=[0 -x(3) x(2); x(3) 0 -x(1); -x(2) x(1) 0];
    E =skew_tij*R_gt;
    err_epipolar(ii) = x2'*E*x1;
    % affine transformation constraints
    At_33(1:2,1:2) = At(:,:,ii);
    At_33(3,3) = 1;
    equationerror = E'*x2 + At_33'*E*x1;
    err_affinetransformation(:,ii) = equationerror(1:2,1); 
end
disp(['The maximum residual of epipolar geometry: ' num2str(max(abs(err_epipolar(:))))])
disp(['The maximum residual of affine transformation constraints: ' num2str(max(abs(err_affinetransformation(:))))])

%% run solver
% 20 solutions
[cay_sols, t_sols, R_sols, cay_sols_all] = solver_depth_mono_2ac(Image1(1:2,:), Image2(1:2,:), At);
cay_sol = find_solution(cay_sols, cay_gt);
t_sol = find_solution_t_singlecamera(t_sols, t_gt);
cay_sol, t_sol, cay_gt, t_gt





