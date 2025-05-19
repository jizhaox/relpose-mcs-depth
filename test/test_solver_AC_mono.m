%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('--------------2AC method & single camera--------------')
n_point = 2;
[Image1, Image2, At, R_gt, cay_gt, t_gt, theta_gt] = generate_2AC_mono_data(n_point);

%% run solver
% 20 solutions
[cay_sols, t_sols, R_sols, cay_sols_all] = solver_depth_mono_2ac(Image1(1:2,:), Image2(1:2,:), At);
cay_sol = find_solution(cay_sols, cay_gt);
t_sol = find_solution_t_mono(t_sols, t_gt);
disp('Cayley: ground truth and estimation');
disp([cay_gt, cay_sol])
disp('Translation: ground truth and estimation');
disp([t_gt, t_sol])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('--------------2AC method & monocular camera with known rotation angle--------------');
disp('--------------depth_mono_2ac_ka--------------');
[cay_sols, t_sols, R_sols, cay_sols_all] = solver_depth_mono_2ac_ka(Image1(1:2,:), Image2(1:2,:), At, theta_gt);
cay_sol = find_solution(cay_sols, cay_gt);
t_sol = find_solution_t_mono(t_sols, t_gt);
disp('Cayley: ground truth and estimation');
disp([cay_gt, cay_sol])
disp('Translation: ground truth and estimation');
disp([t_gt, t_sol])