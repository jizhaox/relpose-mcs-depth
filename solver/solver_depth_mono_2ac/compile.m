if isunix
    % Linux or macOS
   eigen_dir = '/usr/include/eigen3';
elseif ispc
    % Windows
    eigen_dir = 'D:/eigen-3.4.0';
else
    disp('Unknown operating system');
end

tic; mex(['-I"' eigen_dir '"'],'-O','solver_depth_mono_2ac.cpp'); toc
tic; mex(['-I"' eigen_dir '"'],'-O','solver_depth_mono_2ac_ka.cpp'); toc
