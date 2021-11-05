function [data, R_gt, cay_gt, t_gt] = generate_2rpr90_1pt_synthetic()

% Reference:
% [1] Ji Zhao, Laurent Kneip, Yijia He, and Jiayi Ma.
%     Minimal Case Relative Pose Computation using Ray-Point-Ray Features.
%     IEEE Transactions on Pattern Analysis and Machine Intelligence, 
%     42(5): 1176 - 1190, 2020.
% Author: Ji Zhao
% Email: zhaoji84@gmail.com

% Ground truth pose
cay_gt = randn(3,1);
q_gt = [1; cay_gt];
R_gt = quat2rot(normc(q_gt));
t_gt = randn(3,1);
t_gt = t_gt/norm(t_gt);

% setup some 3D points
X = randn(3,3);

% project in first image
x1 = X./X([3;3;3],:);

% create orthogonal directions
D1 = X(:,[1 1]) + 0.1*orth(randn(3,2));
D2 = X(:,[2 2]) + 0.1*orth(randn(3,2));
pd1 = normc(x1(:,[1 1]) - D1./D1([3;3;3],:));
pd2 = normc(x1(:,[2 2]) - D2./D2([3;3;3],:));

d1x = [pd1(:,1) pd2(:,1)];
d1y = [pd1(:,2) pd2(:,2)];

% Change coordinate system to second camera
RX = R_gt*X + repmat(t_gt,[1 size(X,2)]);
RD1 = R_gt*D1 + repmat(t_gt,[1 size(D1,2)]);
RD2 = R_gt*D2 + repmat(t_gt,[1 size(D2,2)]);

% and reproject
x2 = RX./RX([3;3;3],:);
pd1 = normc(x2(:,[1,1]) - RD1./RD1([3;3;3],:));
pd2 = normc(x2(:,[2,2]) - RD2./RD2([3;3;3],:));
d2x = [pd1(:,1) pd2(:,1)];
d2y = [pd1(:,2) pd2(:,2)];

%%
data = [x1(1:2,1), d1x(1:2,1), d1y(1:2,1), x1(1:2,2), d1x(1:2, 2), d1y(1:2, 2), x1(1:2,3), ...
        x2(1:2,1), d2x(1:2,1), d2y(1:2,1), x2(1:2,2), d2x(1:2, 2), d2y(1:2, 2), x2(1:2,3)];
data = data(:);
