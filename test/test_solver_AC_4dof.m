clear
format long

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2AC vertical solver & inter-camera case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('--------------2AC solver with known vertical direction & inter-camera case--------------')
%% generate synthetic data
match_type = 'inter';
[Image1, Image2, At, R_cam, t_cam, R_gt, t_gt, theta_y_gt] = generate_2AC_4dof_data('4DOF', match_type);
%% run solver
[R_sol, t_sol, theta_y_sol] = solver_depth_2ac_4dof(Image1, Image2, At, R_cam, t_cam, match_type);
% display
[~, idx] = min(abs(theta_y_sol - theta_y_gt));
R_sol(:, :, idx), R_gt
t_sol(:, idx), t_gt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2AC vertical solver & intra-camera case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('--------------2AC solver with known vertical direction & intra-camera case--------------')
%% generate synthetic data
match_type = 'intra';
[Image1, Image2, At, R_cam, t_cam, R_gt, t_gt, theta_y_gt] = generate_2AC_4dof_data('4DOF', match_type);
%% run solver
[R_sol, t_sol, theta_y_sol] = solver_depth_2ac_4dof(Image1, Image2, At, R_cam, t_cam, match_type);
% display
[~, idx] = min(abs(theta_y_sol - theta_y_gt));
R_sol(:, :, idx), R_gt
t_sol(:, idx), t_gt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2AC vertical solver & inter-camera case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('--------------2AC solver with known vertical direction & inter-camera case--------------')
%% generate synthetic data
match_type = 'inter';
[Image1, Image2, At, R_cam, t_cam, R_gt, ~, t_gt] = generate_2AC_data(match_type);
[R_imu, theta_y_gt] = generate_imu_measurement(R_gt);

%% run solver
%% transform observations using IMU measurements
[R_sol, t_sol, theta_y_sol] = solver_depth_2ac_4dof_v1(Image1, Image2, At, R_cam, t_cam, R_imu, match_type);
% display
[~, idx] = min(abs(theta_y_sol - theta_y_gt));
R_sol(:, :, idx), R_gt
t_sol(:, idx), t_gt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2AC vertical solver & inter-camera case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('--------------2AC solver with known vertical direction & intra-camera case--------------')
%% generate synthetic data
match_type = 'intra';
[Image1, Image2, At, R_cam, t_cam, R_gt, ~, t_gt] = generate_2AC_data(match_type);
[R_imu, theta_y_gt] = generate_imu_measurement(R_gt);

%% run solver
%% transform observations using IMU measurements
[R_sol, t_sol, theta_y_sol] = solver_depth_2ac_4dof_v1(Image1, Image2, At, R_cam, t_cam, R_imu, match_type);
theta_y_sol,theta_y_gt
% display
[~, idx] = min(abs(theta_y_sol - theta_y_gt));
R_sol(:, :, idx), R_gt
t_sol(:, idx), t_gt

