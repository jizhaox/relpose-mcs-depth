
eigen_dir = '/usr/include/eigen3';

tic; mex(['-I"' eigen_dir '"'],'-O','solver_depth_inter_2ac.cpp'); toc
tic; mex(['-I"' eigen_dir '"'],'-O','solver_depth_intra_2ac.cpp'); toc
