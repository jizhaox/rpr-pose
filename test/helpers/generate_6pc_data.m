function [data0, Tf1tof2, cay_gt] = ...
    generate_6pc_data(match_type)
% Author: Ji Zhao
% Email: zhaoji84@gmail.com
% 01/27/2021
% Reference:
% [1] 

if nargin < 1
    match_type = 1;
end

switch match_type
    case 1
        % inter-camera PC of 2 cameras
        match_info{1} = struct('idx1', 1, 'idx2', 2);
        match_info{2} = struct('idx1', 1, 'idx2', 2);
        match_info{3} = struct('idx1', 1, 'idx2', 2);
        match_info{4} = struct('idx1', 2, 'idx2', 1);
        match_info{5} = struct('idx1', 2, 'idx2', 1);
        match_info{6} = struct('idx1', 2, 'idx2', 1);
    case 3
        % intra-camera PC of 2 cameras
        match_info{1} = struct('idx1', 1, 'idx2', 1);
        match_info{2} = struct('idx1', 1, 'idx2', 1);
        match_info{3} = struct('idx1', 1, 'idx2', 1);
        match_info{4} = struct('idx1', 2, 'idx2', 2);
        match_info{5} = struct('idx1', 2, 'idx2', 2);
        match_info{6} = struct('idx1', 2, 'idx2', 2);
    otherwise
        error('match type error!');
end

%% Part 1: get independent random variables
% generating random extrinsic parameters
n_cam = 2;
avg_cam_distance = 0.5;

cam2f_rotation = cell(n_cam, 1);
cam2f_offset = cell(n_cam, 1);
Rc = cell(n_cam, 1);
tc = cell(n_cam, 1);
for ii = 1:n_cam
    cay = rand(3, p);
    cam2f_rotation{ii} = zp_cayley_rotation(cay, p);
    cam2f_offset{ii} = avg_cam_distance * generateRandomR() * [1.0; 0.0; 0.0];
    
    Rc{ii} = cam2f_rotation{ii}';
    tc{ii} = mod(-Rc{ii}*cam2f_offset{ii}, p);
    % transformation from body reference to perspective camera references
    Hftoc{ii} = [Rc{ii} tc{ii}; 0 0 0 1];
    inv_Hftoc{ii} = mod([Rc{ii}', -Rc{ii}'*tc{ii}; 0, 0, 0, 1], p);
end


% generating random relative pose
[Rf1tof2, quat_gt, u, sin_half_ang, cos_half_ang, tan_half_ang] = zp_quaternion_rotation(p);
cay_gt = mod(u*tan_half_ang, p);
Rf1tof2_check = zp_cayley_rotation(cay_gt, p);
Tf1tof2 = zp_rand(3, p);
% transformation from body reference at time 1 to time 2
Hf1tof2 = [Rf1tof2 Tf1tof2; 0 0 0 1];

%% generating random planes and points
n_point_putative = 10;
point_putative = cell(n_point_putative, 1);
for ii = 1:n_point_putative
    PT = zp_rand(3, p);
    point_putative{ii} = struct('point', PT);
end

%% Part 2: extract point observations
% images at time i
x_i = cell(n_point_putative, n_cam);
for ii = 1:n_point_putative
    PT = point_putative{ii}.point;
    for jj = 1:n_cam
        x_i{ii,jj} = zp_proj_point(mod(Rc{jj}*PT+tc{jj}, p), p);
    end
end
% images at time j
Hftoc_j = cell(n_cam, 1);
Rc_j = cell(n_cam, 1);
tc_j = cell(n_cam, 1);
for ii = 1:n_cam
    tmp = mod(Hftoc{ii}*Hf1tof2, p);
    Hftoc_j{ii} = tmp;
    Rc_j{ii} = tmp(1:3,1:3);
    tc_j{ii} = tmp(1:3,4);
end
x_j = cell(n_point_putative, n_cam);
for ii = 1:n_point_putative
    PT = point_putative{ii}.point;
    for jj = 1:n_cam
        x_j{ii,jj} = zp_proj_point(mod(Rc_j{jj}*PT+tc_j{jj}, p), p);
    end
end

%% Generalized Epipolar Constraint: method 1
epipolar_error = zeros(n_point_putative, 2);
for j = 1:n_point_putative
    for cam_idx = 1:2
        idx1 = match_info{cam_idx}.idx1;
        idx2 = match_info{cam_idx}.idx2;
        % relative pose between time i and time j 
        Hij = mod(Hftoc{idx2}*Hf1tof2*inv_Hftoc{idx1}, p);
        Rij = Hij(1:3,1:3);
        tij = Hij(1:3,4);
        % line observations
        x1 = x_i{j, idx1};
        x2 = x_j{j, idx2};
        E = mod(zp_skew(tij)*Rij, p);
        epipolar_error(j, cam_idx) = mod(x2'*E*x1, p);
    end
end
% check
if ~all(epipolar_error==0)
    error('err_epipolar is wrong!')
end

%% Generalized Epipolar Constraint: method 2
% line observation at time i
Line_i = cell(n_cam, n_point_putative);
for k = 1:n_cam
    for ii = 1:n_point_putative
        [Line_i_tmp, flag] = f_pluckerline(x_i{ii,k}, cam2f_rotation{k}, cam2f_offset{k});
        Line_i{k, ii} = Line_i_tmp;
    end
end
% line observation at time j
Line_j = cell(n_cam, n_point_putative);
for k = 1:n_cam
    for ii = 1:n_point_putative
        [Line_j_tmp, flag] = f_pluckerline(x_j{ii,k}, cam2f_rotation{k}, cam2f_offset{k});
        Line_j{k, ii} = Line_j_tmp;
    end
end

% verification
eEc = zeros(6, 6);
eEc(1:3,1:3) = zp_skew(Tf1tof2)*Rf1tof2;
eEc(1:3,4:6) = Rf1tof2;
eEc(4:6,1:3) = Rf1tof2;
eEc = mod(eEc, p);
eEc_error = zeros(n_point_putative, 2);
for j = 1:n_point_putative
    for cam_idx = 1:2
        idx1 = match_info{cam_idx}.idx1;
        idx2 = match_info{cam_idx}.idx2;
        
        eEc_error(j, cam_idx) = mod(Line_j{idx2,j}'*eEc*Line_i{idx1,j}, p);
    end
end
% check
if ~all(eEc_error==0)
    error('err_epipolar is wrong!')
end

%% define focal length
if (is_calibrated)
    focal_length = 1;
    inv_focal_length = 1;
else
    focal_length = 1;
    flag = false;
    while (focal_length == 1 || flag == false)
        focal_length = zp_rand(1, p);
        [inv_focal_length, flag] = zp_mul_inv(focal_length, p);
    end
end

%% generate data
data0 = [];
num_pc = 6;
for ii = 1:num_pc
    idx_pt = ii;
    cam_idx = ii;

    idx_c1 = match_info{cam_idx}.idx1;
    idx_c2 = match_info{cam_idx}.idx2;
    d = [cam2f_rotation{idx_c1}(:); cam2f_offset{idx_c1}(:);
         cam2f_rotation{idx_c2}(:); cam2f_offset{idx_c2}(:);
         mod(x_i{idx_pt,idx_c1}(1:2)*focal_length, p);
         mod(x_j{idx_pt,idx_c2}(1:2)*focal_length, p)];
    data0 = [data0; d(:)];
end

if is_known_angle
    data0 = [data0; tan_half_ang];
end