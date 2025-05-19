function [err_epipolar, err_affinetransformation] = check_constraint(Image1, Image2, At, R_cam, t_cam, R_gt, t_gt, match_info)

%% check data
% For noise-free data, the residual of epipolar constraint and affine transformation constraints should be small (Eqs.(7) and (8) in the paper)
err_epipolar = zeros(2, 1);
err_affinetransformation = zeros(2, 2);
Tf_gt = [R_gt t_gt; 0 0 0 1];
Tf_cam{1} = [R_cam(:,:,1), t_cam(:,1); 0, 0, 0, 1];
Tf_cam{2} = [R_cam(:,:,2), t_cam(:,2); 0, 0, 0, 1];
for ii = 1:2
    idx1 = match_info{ii}.idx1;
    idx2 = match_info{ii}.idx2;
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