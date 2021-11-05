function [L1, L2, v1c1_v2c2, v1c2_v2c1, t, R ] = create2D2DExperiment_inter( pt_number, noise )

cam_number = 2;
v1c1_v2c2 = zeros(6, pt_number); % correspondence between view1_cam1 and view2_cam2
v1c2_v2c1 = zeros(6, pt_number); % correspondence between view1_cam2 and view2_cam1

%% generate the camera system

avg_cam_distance = 0.5;
cam_offsets = zeros(3,cam_number);
%cam_rotations = zeros(3,cam_number*3);

for i=1:cam_number
    cam_offsets(:,i) = avg_cam_distance * generateRandomR() * [1.0; 0.0; 0.0];
    %cam_rotations(:,(i-1)*3+1:(i-1)*3+3) = generateRandomR();
end


%% generate random view-points

max_parallax = 2.0;
max_rotation = 0.5;

position1 = zeros(3,1);
rotation1 = eye(3);

position2 = max_parallax * 2.0 * (rand(3,1) - repmat(0.5,3,1));
rotation2 = generateBoundedR(max_rotation);

%% Generate random point-cloud

minDepth = 4.0;
maxDepth = 8.0;

normalizedPoints = 2.0*(rand(3,pt_number)-repmat(0.5,3,pt_number));
norms = sqrt(sum(normalizedPoints.*normalizedPoints));
directions = normalizedPoints./repmat(norms,3,1);
points = (maxDepth-minDepth) * normalizedPoints + minDepth * directions;

%% Now create the correspondences by looping through the cameras

focal_length = 800.0;

L1 = zeros(6,pt_number);
L2 = zeros(6,pt_number);
cam_correspondence = 1;
cam_correspondences = zeros(1,pt_number);

for i=1:pt_number
    
    cam_offset = cam_offsets(:,cam_correspondence);
    %cam_rotation = cam_rotations(:,(cam_correspondence-1)*3+1:(cam_correspondence-1)*3+3);
    
    body_point1 = rotation1' * (points(:,i)-position1);
    body_point2 = rotation2' * (points(:,i)-position2);
    
    % we actually omit the can rotation here by unrotating the bearing
    % vectors already
    bearingVector1 = body_point1 - cam_offset;
    bearingVector2 = body_point2 - cam_offset;
    bearingVector1_norm = norm(bearingVector1);
    bearingVector2_norm = norm(bearingVector2);
    bearingVector1 = bearingVector1/bearingVector1_norm;
    bearingVector2 = bearingVector2/bearingVector2_norm;
    
    % add noise to the bearing vectors here
    bearingVector1_noisy = addNoise(bearingVector1,focal_length,noise);
    bearingVector2_noisy = addNoise(bearingVector2,focal_length,noise);
    
    % store the normalized bearing vectors along with the cameras they are
    % being seen (we create correspondences that always originate from the
    % same camera, you can change this if you want)
    bearingVector1_norm = norm(bearingVector1_noisy);
    bearingVector2_norm = norm(bearingVector2_noisy);
    
    % plucker line representation
    tmp1 = bearingVector1_noisy./bearingVector1_norm;
    L1(:,i) = [tmp1; cross(cam_offset,tmp1)];
    tmp2 = bearingVector2_noisy./bearingVector2_norm;
    L2(:,i) = [tmp2; cross(cam_offset,tmp2)];

    % change the camera correspondence
    cam_correspondences(1,i) = cam_correspondence;
    cam_correspondence = cam_correspondence + 1;
    if cam_correspondence > cam_number
        cam_correspondence = 1;
    end
end

%% compute relative translation and rotation

R = rotation1' * rotation2;
t = rotation1' * (position2 - position1);


