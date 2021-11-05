%% test solvers for 2RPR90+1pt method

% Reference:
% [1] Ji Zhao, Laurent Kneip, Yijia He, and Jiayi Ma.
%     Minimal Case Relative Pose Computation using Ray-Point-Ray Features.
%     IEEE Transactions on Pattern Analysis and Machine Intelligence, 
%     42(5): 1176 - 1190, 2020.
% Author: Ji Zhao
% Email: zhaoji84@gmail.com

%% Interface of the solver
% [R_sols, t_sols, cay_sols] = solver_2rpr90_1pt(data);
% INPUT
% data: RPR90 and point observations of two views
%     size is 28*1
%     order is: RPR90 #1, RPR #2, point #1 at view 1, 
%               RPR90 #1, RPR #2, point #1 at view 2.
%     Each RPR90 structure includes normalized image coordinates and 
%     direction vectors of two rays on the image plane.
% OUTPUT
% R_sols: real solutions for rotation using rotation matrix
%     size is 3*3*N, where N is the number of real solutions
% t_sols: real solutions for translation
%     size is 3*N, where N is the number of real solutions
% IMPORTANT NOTE: Due to the scale ambiguity of translation, we normalize
%     the resulted translation vectors. Still, there remains sign ambiguity
%     of the translation vectors. When using hypothesis-and-test framework
%     to find the best rotation-translation pair, you should try
%     both of the transltion vectors and their opposite directions.
% cay_sols: real solutions for rotation using Cayley
%     size is 3*N, where N is the number of real solutions

clear;
disp('================== 2RPR90 + 1pt method ==================');
%% prepare data
[data, R_gt, cay_gt, t_gt] = generate_2rpr90_1pt_synthetic();

%% run the solver
tic, [R_sols, t_sols] = solver_2rpr90_1pt(data); toc
R_sols, R_gt
t_sols, t_gt

