function [data, R_gt, cay_gt, t_gt] = generate_1rpr90_3pt_synthetic()

% Reference:
% [1] Ji Zhao, Laurent Kneip, Yijia He, and Jiayi Ma.
%     Minimal Case Relative Pose Computation using Ray-Point-Ray Features.
%     IEEE Transactions on Pattern Analysis and Machine Intelligence, 
%     42(5): 1176 - 1190, 2020.
% Author: Ji Zhao
% Email: zhaoji84@gmail.com

c = inf;
thresh = cos(deg2rad(45)); % threshold to reject near-degenerate cases

while(c>thresh)
    % Ground truth pose
    q_gt = rand(4, 1)+0.1;
    cay_gt = q_gt(2:4)/q_gt(1);
    R_gt = quat2rot(normc(q_gt));
    t_gt = 1+rand(3,1); % length of translation is [1, 2]

    % setup some 3D points
    X = (1+rand(3,4))*5; % point depth is [5, 10]

    % project in first image
    x1 = X./X([3;3;3],:);

    % create orthogonal directions
    D1 = X(:,[1 1]) + 0.1*orth(randn(3,2));
    pd1 = normc(x1(:,[1 1]) - D1./D1([3;3;3],:));

    d1x = pd1(:,1); d1x = normc(d1x);
    d1y = pd1(:,2); d1y = normc(d1y);

    % Change coordinate system to second camera
    RX = R_gt*X + repmat(t_gt,[1 size(X,2)]);
    RD1 = R_gt*D1 + repmat(t_gt,[1 size(D1,2)]);

    % and reproject
    x2 = RX./RX([3;3;3],:);
    pd1 = normc(x2(:,[1,1]) - RD1./RD1([3;3;3],:));
    d2x = pd1(:,1); d2x = normc(d2x);
    d2y = pd1(:,2); d2y = normc(d2y);

    c1 = abs(d1x'*d1y/(norm(d1x)*norm(d1y)));
    c2 = abs(d2x'*d2y/(norm(d2x)*norm(d2y)));
    c = max([c1; c2]);
end
%%
data = [x1(1:2,1), d1x(1:2,1), d1y(1:2,1), x1(1:2,2), x1(1:2,3), x1(1:2,4), ...
        x2(1:2,1), d2x(1:2,1), d2y(1:2,1), x2(1:2,2), x2(1:2,3), x2(1:2,4)];
data = data(:);

t_gt = t_gt/norm(t_gt);