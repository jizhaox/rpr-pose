
%eigen_dir = 'C:\jizhao\project\Eigen'; 
eigen_dir = '/usr/include/eigen3';

tic; mex(['-I"' eigen_dir '"'],'-O','solver_1rpr90_1pt.cpp'); toc
tic; mex(['-I"' eigen_dir '"'],'-O','solver_1rpr90_3pt.cpp'); toc
tic; mex(['-I"' eigen_dir '"'],'-O','solver_2rpr90_1pt.cpp'); toc
tic; mex(['-I"' eigen_dir '"'],'-O','solver_fivept.cpp'); toc
