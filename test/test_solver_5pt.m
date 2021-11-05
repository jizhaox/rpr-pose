%% test solvers for five-point method

% Reference:
% [1] Ji Zhao, Laurent Kneip, Yijia He, and Jiayi Ma.
%     Minimal Case Relative Pose Computation using Ray-Point-Ray Features.
%     IEEE Transactions on Pattern Analysis and Machine Intelligence, 
%     42(5): 1176 - 1190, 2020.
% Author: Ji Zhao
% Email: zhaoji84@gmail.com

%% Interface of the solver
% [R_sols, t_sols, cay_sols] = solver_fivept(data);
% INPUT
% data: normalized image coordinates of two views
%     size is 20*1
%     order is: points 1,2,3,4,5 at view 1, points 1,2,3,4,5 at view 2 
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
disp('================== 5 point method ==================');
%% generate synthetic data by openGV
cam_number = 1;
pt_number = 50;
outlier_fraction = 0.0;
noise = 0.0;

[P1,P2,t_gt,R_gt] = create2D2DExperiment(pt_number,cam_number,noise,outlier_fraction);
t_gt = t_gt/norm(t_gt(:));

%% run the solver
P1_img = P1(1:2,:)./repmat(P1(3,:), [2 1]);
P2_img = P2(1:2, :)./repmat(P2(3,:), [2 1]);
pt1 = P1_img(:, 1:5);
pt2 = P2_img(:, 1:5);
data = [pt2(:); pt1(:)];
tic, [R_sols, t_sols] = solver_fivept(data); toc
R_sols, R_gt
t_sols, t_gt

