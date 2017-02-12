function [C, R] = LinearPnP(X, x, K)
%% LinearPnP
% Getting pose from 2D-3D correspondences
% Inputs:
%     X - size (N x 3) matrix of 3D points
%     x - size (N x 2) matrix of 2D points whose rows correspond with X
%     K - size (3 x 3) camera calibration (intrinsics) matrix
% Outputs:
%     C - size (3 x 1) pose transation
%     R - size (3 x 1) pose rotation
%
% IMPORTANT NOTE: While theoretically you can use the x directly when solving
% for the P = [R t] matrix then use the K matrix to correct the error, this is
% more numeically unstable, and thus it is better to calibrate the x values
% before the computation of P then extract R and t directly

% Philippe Weingerner January 2017

[N, ~] = size(X);

A = [];

for i = 1:N
    Xi = [X(i, :), 1]; % 1x4
    Ai = zeros(3, 12);
    Ai(1, 1:4) = Xi;
    Ai(2, 5:8) = Xi;
    Ai(3, 9:12) = Xi;
    xc = inv(K) * [x(i, 1); x(i, 2); 1];
    Ai = Vec2Skew(xc) * Ai;
    A = [A; Ai];
end

[U, D, V] = svd(A);
P = reshape(V(:,end), [4,3])';

%R = K' * P(:, 1:3);
R = P(:, 1:3);

[U, D, V] = svd(R);
R = U * V'; % enforcing orthogonal constraint of rotation matrix

%t = K' * P(:, 4) ./ D(1,1); % translation vector
t = P(:, 4) ./ D(1,1); % translation vector

C = -R' * t; % Camera Optical Center Coordinates

if det(R) < 0 % == -1
    R = -R;
    C = -C;
end





