function [err_epipolar, err_affinetransformation] = check_constraint_mono(Image1, Image2, At, R_gt, t_gt)

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