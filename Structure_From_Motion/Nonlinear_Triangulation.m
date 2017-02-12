function X = Nonlinear_Triangulation(K, C1, R1, C2, R2, C3, R3, x1, x2, x3, X0)
%% Nonlinear_Triangulation
% Refining the poses of the cameras to get a better estimate of the points
% 3D position
% Inputs: 
%     K  - size (3 x 3) camera calibration (intrinsics) matrix
%     Ci - size (3 x 1)
%     Ri - size (3 x 3)
%     xi - size (N x 2)
% Outputs: 
%     X - size (N x 3) matrix of refined point 3D locations 

% Philippe Weingerner January 2017

[N, ~] = size(X0);
X = zeros(N, 3);

% Per 3D point processing:
% unknwonws (X0 to be fine tuned): 3 (3D coordinates)
% outputs: 6 (2D coordinates for 3 different views of above 3D point) 
% f: R3 -> R6

% Jacobian: 6x3

for i = 1:N
    X(i, :) = Single_Point_Nonlinear_Triangulation(K, C1, R1, C2, R2, C3, R3, x1(i,:), x2(i,:), x3(i,:), X0(i, :));
end

end

function fX = f(K, C1, R1, C2, R2, C3, R3, X)
    xp1 = K * R1 * (X - C1);
    xp2 = K * R2 * (X - C2);
    xp3 = K * R3 * (X - C3);
    xp1 = xp1 ./ xp1(3);
    xp2 = xp2 ./ xp2(3);
    xp3 = xp3 ./ xp3(3);
    fX = [xp1(1:2); xp2(1:2); xp3(1:2)]; % 6x1
end

function X = Single_Point_Nonlinear_Triangulation(K, C1, R1, C2, R2, C3, R3, x1, x2, x3, X0)
    b = [ x1 x2 x3 ]'; % 6x1
    %Xi = X0(i, :)'; % 3x1
    Xi = X0';
    X  = X0';
    fX = f(K, C1, R1, C2, R2, C3, R3, Xi);
    error = norm(fX-b);
    
    %while (error_decrease > 0)
    for cnt=1:1
        J1 = Jacobian_Triangulation(C1, R1, K, Xi); % 2x3
        J2 = Jacobian_Triangulation(C2, R2, K, Xi); % 2x3
        J3 = Jacobian_Triangulation(C3, R3, K, Xi); % 2x3
        J = [J1; J2; J3]; % 6x3
        %if rcond(J'*J) < 1e-4
        %    break;
        %end
        dX = (J'*J) \ (J' * (b -fX)); % 3x1 : inv(3x6 6x3) * 3x6 * (6x1)
        Xi = Xi + dX;
        fX = f(K, C1, R1, C2, R2, C3, R3, Xi);
        error = norm(fX-b);
    end
    X = Xi;
end

function J = Jacobian_Triangulation(C, R, K, X)
    KR = K * R; % 3x3
    du = KR(1, :); % 1x3
    dv = KR(2, :); % 1x3
    dw = KR(3, :); % 1x3
    xp = KR * (X - C); % 3x1
    u = xp(1); % 1x1
    v = xp(2); % 1x1
    w = xp(3); % 1x1
    J = [(w.*du-u.*dw)./w^2; (w.*dv-v.*dw)./w^2]; % 2x3
end
