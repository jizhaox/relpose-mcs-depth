function [R_imu, theta_y_gt] = generate_imu_measurement(R)

% Convert rotation matrix to Euler angles, and assume only yaw is unknown
eul = rotm2eul(R, 'ZYX');
theta_z = eul(1);
theta_y = eul(2);
theta_x = eul(3);
R_left = [cos(theta_z), -sin(theta_z), 0;
      sin(theta_z), cos(theta_z), 0;
      0, 0, 1];
Ry = [cos(theta_y), 0, sin(theta_y);
      0, 1, 0;
     -sin(theta_y), 0, cos(theta_y)];
R_right = [1, 0, 0;
      0, cos(theta_x), -sin(theta_x);
      0, sin(theta_x), cos(theta_x)];
R_imu = cat(3, R_left, R_right);

% There is a minus sign to align with our definition of Ry
theta_y_gt = -theta_y;

% check
%R_err = R_left * Ry * R_right - R
