function F = EstimateFundamentalMatrix(x1, x2)
%% EstimateFundamentalMatrix
% Estimate the fundamental matrix from two image point correspondences 
% Inputs:
%     x1 - size (N x 2) matrix of points in image 1
%     x2 - size (N x 2) matrix of points in image 2, each row corresponding
%       to x1
% Output:
%    F - size (3 x 3) fundamental matrix with rank 2

% Philippe Weingerner January 2017

[N, ~] = size(x1);
A = zeros(N, 3*3);

for i = 1:N
    % !!! there is a typo in the course slides x1i <-> x2i !!!
    x1i = [x1(i, 1) x1(i, 2) 1]; % 1x3
    %x2i = [x2(i, 1) x2(i, 2) 1];
    A(i, :) = [x2(i,1).*x1i, x2(i,2).*x1i, x1i]; % 1x9 or kron(x1i, x2i)
    %A(i, :) = kron(x2i, x1i);
end;

[~, ~, V] = svd(A);
F = reshape(V(:,end), [3,3])';

[U, D, V] = svd(F);
D(3,3) = 0; % make sure F ends up being rank 2
F = U * D * V';

%norm(F) should be 1

