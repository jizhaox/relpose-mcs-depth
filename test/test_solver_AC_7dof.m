clear;
format long

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('--------------3AC method with unknow scale & inter-camera case--------------')
%% generate synthetic data
num_ac = 3;
[Image1, Image2, At, R_cam, t_cam, R_gt, cay_gt, t_gt, ~, match_info] = generate_2AC_data('inter', num_ac);
%% simulate unknown scale
s_gt = 0.5 + rand(1);
R_cam_first = zeros(3, 3, num_ac);
R_cam_second = zeros(3, 3, num_ac);
t_cam_first = zeros(3, num_ac);
t_cam_second = zeros(3, num_ac);
for ii = 1:num_ac
    idx1 = match_info{ii}.idx1;
    idx2 = match_info{ii}.idx2;
    R_cam_first(:, :, ii) = R_cam(:, :, idx1);
    R_cam_second(:, :, ii) = R_cam(:, :, idx2);
    t_cam_first(:, ii) = t_cam(:, idx1);
    % introduce the unknown scale
    t_cam_second(:, ii) = t_cam(:, idx2) / s_gt;
end

%% run solver
[cay_sols, t_sols, R_sols, s_sols] = solver_depth_3ac(Image1, Image2, At, R_cam_first, R_cam_second, t_cam_first, t_cam_second);
[cay_sol, idx] = find_solution(cay_sols, cay_gt);
t_sol = t_sols(:, idx);
s_sol = s_sols(idx);
disp('Cayley: ground truth and estimation');
disp([cay_gt, cay_sol])
disp('Translation: ground truth and estimation');
disp([t_gt, t_sol])
disp('scale: ground truth and estimation');
disp([s_gt, s_sol])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('--------------3AC method with unknow scale & intra-camera case--------------')
%% generate synthetic data
num_ac = 3;
[Image1, Image2, At, R_cam, t_cam, R_gt, cay_gt, t_gt, ~, match_info] = generate_2AC_data('intra', num_ac);
%% simulate unknown scale
s_gt = 0.5 + rand(1);
R_cam_first = zeros(3, 3, num_ac);
R_cam_second = zeros(3, 3, num_ac);
t_cam_first = zeros(3, num_ac);
t_cam_second = zeros(3, num_ac);
for ii = 1:num_ac
    idx1 = match_info{ii}.idx1;
    idx2 = match_info{ii}.idx2;
    R_cam_first(:, :, ii) = R_cam(:, :, idx1);
    R_cam_second(:, :, ii) = R_cam(:, :, idx2);
    t_cam_first(:, ii) = t_cam(:, idx1);
    % introduce the unknown scale
    t_cam_second(:, ii) = t_cam(:, idx2) / s_gt;
end
%% run solver
[cay_sols, t_sols, R_sols, s_sols] = solver_depth_3ac(Image1, Image2, At, R_cam_first, R_cam_second, t_cam_first, t_cam_second);
[cay_sol, idx] = find_solution(cay_sols, cay_gt);
t_sol = t_sols(:, idx);
s_sol = s_sols(idx);
disp('Cayley: ground truth and estimation');
disp([cay_gt, cay_sol])
disp('Translation: ground truth and estimation');
disp([t_gt, t_sol])
disp('scale: ground truth and estimation');
disp([s_gt, s_sol])
