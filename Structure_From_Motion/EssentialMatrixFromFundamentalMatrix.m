function E = EssentialMatrixFromFundamentalMatrix(F,K)
%% EssentialMatrixFromFundamentalMatrix
% Use the camera calibration matrix to esimate the Essential matrix
% Inputs:
%     K - size (3 x 3) camera calibration (intrinsics) matrix
%     F - size (3 x 3) fundamental matrix from EstimateFundamentalMatrix
% Outputs:
%     E - size (3 x 3) Essential matrix with singular values (1,1,0)

% Philippe Weingerner January 2017

E = K' * F * K;
[U, D, V] = svd(E);
E = U * diag([1,1,0]) * V'; % make sure E is rank 2 with norm 1
