function [Image1, Image2, At, R_extrinsic_para, T_extrinsic_para, R_gt, cay_gt, t_gt, theta_gt, match_info] = generate_2AC_data(match_type)

%% Two types of affine correspondence
%% inter-cam affine correspondences contain
% (1) 1 affine correspondence: (cam1 at time i) <--> (cam2 at time j)
% (2) 1 affine correspondence: (cam2 at time i) <--> (cam1 at time j)
%% intra-cam affine correspondences contain
% (1) 1 affine correspondence: (cam1 at time i) <--> (cam1 at time j)
% (2) 1 affine correspondence: (cam2 at time i) <--> (cam2 at time j)

match_info = match_type_gcam_ac(match_type);
n_cam = 2;
n_point = 2;

%% define random extrinsic parameters
cam_body_rotation = cell(n_cam, 1);
cam_body_offset = cell(n_cam, 1);
R_cam = cell(n_cam, 1);
t_cam = cell(n_cam, 1);
T_body_cam = cell(n_cam, 1);
for ii = 1:n_cam
    cay = rand(3, 1);
    cam_body_rotation{ii} = cayley_rotation(cay);
    cam_body_offset{ii} = rand(3, 1);
    
    R_cam{ii} = cam_body_rotation{ii}';
    t_cam{ii} = -R_cam{ii}*cam_body_offset{ii};
    % transformation from body reference to perspective camera references
    T_body_cam{ii} = [R_cam{ii} t_cam{ii}; 0 0 0 1];
end

%% define relative pose
cay = rand(3, 1);
R_gt = cayley_rotation(cay);
q = rotm2quat(R_gt);
cay_gt = q(2:4)/q(1);
cay_gt = cay_gt(:);
theta_gt = acos(2*q(1)^2-1);

t_gt = rand(3, 1);

% transformation from body reference at time i to time j
T_gt = [R_gt t_gt; 0 0 0 1];


%% generating affine correspondences
% generating random scene points
[PT, Distance, Nvector] = generate_3Dscenepoints(n_point);
points_all = cell(n_point, 1);
for ii = 1:n_point
    points_all{ii} = struct('point', PT(:,:,ii));
end

%% extract point observations
% images at time i
x_i = cell(n_point, n_cam);
for ii = 1:n_point
    PT = points_all{ii}.point;
    for jj = 1:n_cam
        tmp = R_cam{jj}*PT+t_cam{jj};
        x_i{ii,jj} = tmp./tmp(3,:);
    end
end
% images at time j
Rc_j = cell(n_cam, 1);
tc_j = cell(n_cam, 1);
for ii = 1:n_cam
    tmp = T_body_cam{ii}*T_gt;
    Rc_j{ii} = tmp(1:3,1:3);
    tc_j{ii} = tmp(1:3,4);
end
x_j = cell(n_point, n_cam);
for ii = 1:n_point
    PT = points_all{ii}.point;
    for jj = 1:n_cam
        tmp = Rc_j{jj}*PT+tc_j{jj};
        x_j{ii,jj} = tmp./tmp(3,:);
    end
end

%% construct observations
R_extrinsic_para = cat(3, cam_body_rotation{1}, cam_body_rotation{2});
T_extrinsic_para = [cam_body_offset{1}, cam_body_offset{2}];
Image1 = zeros(3, 2);
Image2 = zeros(3, 2);
At = zeros(2, 2, 2);

for ii = 1:n_point
    idx1 = match_info{ii}.idx1;
    idx2 = match_info{ii}.idx2;

    Image1(:,ii) = x_i{ii,idx1}(:,1);
    Image2(:,ii) = x_j{ii,idx2}(:,1);
    H_c1i_c2j = DLT_homography(x_i{ii,idx1}(:,2:end)',x_j{ii,idx2}(:,2:end)');
    At(:,:,ii) = get_affinetransformation(H_c1i_c2j, Image1(:,ii), Image2(:,ii));
end

check_constraint(Image1, Image2, At, R_extrinsic_para, T_extrinsic_para, R_gt, t_gt, match_info);
