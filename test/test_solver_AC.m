clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('--------------2AC method & inter-camera case--------------')
%% generate synthetic data
[Image1, Image2, At, R_cam, t_cam, R_gt, cay_gt, t_gt, theta_gt] = generate_2AC_data('inter');
    
%% run solver
% set the 6th parameter as 0 (default), 48 solutions
[cay_sols, t_sols, R_sols, cay_sols_all] = solver_depth_inter_2ac(Image1, Image2, At, R_cam, t_cam, 0);
cay_sol = find_solution(cay_sols, cay_gt);
t_sol = find_solution(t_sols, t_gt);
disp('Cayley: ground truth and estimation');
disp([cay_gt, cay_sol])
disp('Translation: ground truth and estimation');
disp([t_gt, t_sol])

% set the 6th parameter as 1, 56 solutions.
% according to our evaluation, 56-solution configuration has better numerical stability
[cay_sols, t_sols, R_sols, cay_sols_all] = solver_depth_inter_2ac(Image1, Image2, At, R_cam, t_cam, 1);
cay_sol = find_solution(cay_sols, cay_gt);
t_sol = find_solution(t_sols, t_gt);
disp('Cayley: ground truth and estimation');
disp([cay_gt, cay_sol])
disp('Translation: ground truth and estimation');
disp([t_gt, t_sol])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('--------------2AC method & intra-camera case--------------')
%% generate synthetic data
[Image1, Image2, At, R_cam, t_cam, R_gt, cay_gt, t_gt, theta_gt] = generate_2AC_data('intra');

%% run solver
% 48 solutions
[cay_sols, t_sols, R_sols, cay_sols_all] = solver_depth_intra_2ac(Image1, Image2, At, R_cam, t_cam);
cay_sol = find_solution(cay_sols, cay_gt);
t_sol = find_solution(t_sols, t_gt);
disp('Cayley: ground truth and estimation');
disp([cay_gt, cay_sol])
disp('Translation: ground truth and estimation');
disp([t_gt, t_sol])

% add coordinate transformation
% according to our evaluation, it has better numerical stability for some cases
[cay_sols, t_sols, R_sols] = solver_depth_intra_2ac_enhance(Image1, Image2, At, R_cam, t_cam);
cay_sol = find_solution(cay_sols, cay_gt);
t_sol = find_solution(t_sols, t_gt);
disp('Cayley: ground truth and estimation');
disp([cay_gt, cay_sol])
disp('Translation: ground truth and estimation');
disp([t_gt, t_sol])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('--------------2AC method & inter-camera case with known rotation angle--------------');
%% generate synthetic data
[Image1, Image2, At, R_cam, t_cam, R_gt, cay_gt, t_gt, theta_gt] = generate_2AC_data('inter');

%% run solver
% default, 36 solutions
[cay_sols, t_sols, R_sols, cay_sols_all] = solver_depth_inter_2ac_ka(Image1, Image2, At, R_cam, t_cam, theta_gt, 0);
cay_sol = find_solution(cay_sols, cay_gt);
t_sol = find_solution(t_sols, t_gt);
disp('Cayley: ground truth and estimation');
disp([cay_gt, cay_sol])
disp('Translation: ground truth and estimation');
disp([t_gt, t_sol])

% set the 7th parameter as non-zero, 42 solutions
[cay_sols, t_sols, R_sols, cay_sols_all] = solver_depth_inter_2ac_ka(Image1, Image2, At, R_cam, t_cam, theta_gt, 1);
cay_sol = find_solution(cay_sols, cay_gt);
t_sol = find_solution(t_sols, t_gt);
disp('Cayley: ground truth and estimation');
disp([cay_gt, cay_sol])
disp('Translation: ground truth and estimation');
disp([t_gt, t_sol])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('--------------2AC method & intra-camera case with known rotation angle--------------');
%% generate synthetic data
[Image1, Image2, At, R_cam, t_cam, R_gt, cay_gt, t_gt, theta_gt] = generate_2AC_data('intra');
%% run solver
% default, 36 solutions
[cay_sols, t_sols, R_sols, cay_sols_all] = solver_depth_intra_2ac_ka(Image1, Image2, At, R_cam, t_cam, theta_gt, 0);
cay_sol = find_solution(cay_sols, cay_gt);
t_sol = find_solution(t_sols, t_gt);
disp('Cayley: ground truth and estimation');
disp([cay_gt, cay_sol])
disp('Translation: ground truth and estimation');
disp([t_gt, t_sol])

% set the 7th parameter as non-zero, 44 solutions
[cay_sols, t_sols, R_sols, cay_sols_all] = solver_depth_intra_2ac_ka(Image1, Image2, At, R_cam, t_cam, theta_gt, 1);
cay_sol = find_solution(cay_sols, cay_gt);
t_sol = find_solution(t_sols, t_gt);
disp('Cayley: ground truth and estimation');
disp([cay_gt, cay_sol])
disp('Translation: ground truth and estimation');
disp([t_gt, t_sol])
