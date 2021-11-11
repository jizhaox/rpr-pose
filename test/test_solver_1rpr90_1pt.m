% test solvers for monocular cameras with ray-point-ray features and 
% known rotation axis

% Reference:
% [1] Ji Zhao, Laurent Kneip, Yijia He, and Jiayi Ma.
%     Minimal Case Relative Pose Computation using Ray-Point-Ray Features.
%     IEEE Transactions on Pattern Analysis and Machine Intelligence, 
%     42(5): 1176 - 1190, 2020.
% Author: Ji Zhao
% Email: zhaoji84@gmail.com
% https://sites.google.com/site/drjizhao/

%% Interface of the solver
% [R_sols, t_sols, q_sols] = solver_1rpr90_1pt(data, rotation_axis);
% INPUT
% data: RPR90 and point observations of two views
%     size is 16*1
%     order is: RPR90 #1, point #1 at view 1, 
%               RPR90 #1, point #1 at view 2.
%     Each RPR90 structure includes normalized image coordinates and 
%     direction vectors of two rays on the image plane.
% rotation_axis: rotation axis, support 'x', 'y', and 'z'
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
% q_sols: real solutions for tan(rotation_angle/2)
%     size is 1*N, where N is the number of real solutions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('================== 1RPR90 + 1pt with known rotation axis ==================');
rotation_axis = 'y';
[data, R_gt, theta, t_gt] = generate_1rpr90_1pt_synthetic(rotation_axis);

[R_sols, t_sols, q_sols] = solver_1rpr90_1pt(data, rotation_axis);
R_sols, R_gt
t_sols, t_gt
