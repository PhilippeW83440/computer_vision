function X = LinearTriangulation(K, C1, R1, C2, R2, x1, x2)
%% LinearTriangulation
% Find 3D positions of the point correspondences using the relative
% position of one camera from another
% Inputs:
%     C1 - size (3 x 1) translation of the first camera pose
%     R1 - size (3 x 1) rotation of the first camera pose
%     C2 - size (3 x 1) translation of the second camera
%     R2 - size (3 x 1) rotation of the second camera pose
%     x1 - size (N x 2) matrix of points in image 1
%     x2 - size (N x 2) matrix of points in image 2, each row corresponding
%       to x1
% Outputs: 
%     X - size (N x 3) matrix whos rows represent the 3D triangulated
%       points

% Philippe Weingerner January 2017

%vec2skew = @(v) [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];

% Projection matrices
P1 = K * R1 * [eye(3) -C1];
P2 = K * R2 * [eye(3) -C2];

[N, ~] = size(x1);
X = zeros(N, 3);

for i=1:N
    x1i = [x1(i, 1) x1(i, 2) 1]; % 1x3
    x2i = [x2(i, 1) x2(i, 2) 1];
    % we will solve A * X = 0 in a Linear Least Square way
    A = [Vec2Skew(x1i) * P1; Vec2Skew(x2i) * P2];
    [~, ~, V] = svd(A);
    Xi = V(:, end) / V(end, end);
    X(i, :) = Xi(1:3);
end




