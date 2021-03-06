function [ H ] = est_homography(video_pts, logo_pts)
% est_homography estimates the homography to transform each of the
% video_pts into the logo_pts
% Inputs:
%     video_pts: a 4x2 matrix of corner points in the video
%     logo_pts: a 4x2 matrix of logo points that correspond to video_pts
% Outputs:
%     H: a 3x3 homography matrix such that logo_pts ~ H*video_pts

% Philippe WEINGERTNER december 2016

%H = [];
%H = zeros(3, 3);

[N, dim]= size(video_pts);
A = zeros(2*N, 3*3);

for i=1:N
    x = video_pts(i, :);
    xp= logo_pts(i, :);
    A(2*i-1, :) = [-x(1), -x(2), -1, 0, 0, 0, x(1)*xp(1), x(2)*xp(1), xp(1)];
    A(2*i,   :) = [0, 0, 0, -x(1), -x(2), -1, x(1)*xp(2), x(2)*xp(2), xp(2)];
end

[U, S, V] = svd(A);
H = reshape(V(:,end), [3,3])';

end

