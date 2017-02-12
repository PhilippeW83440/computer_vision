
function [t, R] = ar_test(H)
%% ar_cube
% Estimate your position and orientation with respect to a set of 4 points on the ground
% Inputs:
%    H - the computed homography from the corners in the image
%    render_points - size (N x 3) matrix of world points to project
%    K - size (3 x 3) calibration matrix for the camera
% Outputs: 
%    proj_points - size (N x 2) matrix of the projected points in pixel
%      coordinates
%    t - size (3 x 1) vector of the translation of the transformation
%    R - size (3 x 3) matrix of the rotation of the transformation


% Philippe Weingertner January 2017: Extract the pose from the homography
if H(3,3) < 0
    H = -H;
end
Rp = [ H(:,1), H(:,2), cross(H(:,1), H(:,2)) ];
[U, S, V] = svd(Rp);
R = U * diag([1, 1, det(U*V')]) * V'; 
t = H(:,3) / norm(H(:,1));
%t = H(:,3) / (0.5 * (norm(H(:,1)) + norm(H(:,2))));


end


