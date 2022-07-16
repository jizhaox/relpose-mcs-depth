
eigen_dir = '/usr/include/eigen3';

tic; mex(['-I"' eigen_dir '"'],'-O','solver_depth_mono_2ac.cpp'); toc

