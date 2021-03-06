function [ corners ] = track_corners(images, img_pts_init)
%TRACK_CORNERS 
% This function tracks the corners in the image sequence and visualizes a
% virtual box projected into the image
% Inputs:
%     images - size (N x 1) cell containing the sequence of images to track
%     img_pts_init - size (4 x 2) matrix containing points to initialize
%       the tracker
% Outputs:
%     corners - size (4 x 2 x N) array of where the corners are tracked to

corners = zeros(4,2,size(images,1));

% Philippe Weingertner January 2017
%%%% INITIALIZATION CODE FOR TRACKER HERE %%%%

% Create a point tracker 
pointTracker = vision.PointTracker('MaxBidirectionalError', 2);
% Initialize the tracker
initialize(pointTracker, img_pts_init, images{1});

img_pts = img_pts_init; % img_pts is where you will store the tracked points
corners(:,:,1) = img_pts;

% Iterate through the rest of the images
for i = 2:size(images,1)
    %%%% CODE FOR TRACKING HERE %%%%
    
    % Track the points. Note that some points may be lost.
    [img_pts, isFound] = step(pointTracker, images{i});
    %fprintf('%d\n', isFound);
    if sum(isFound) < 4
        fprintf('Not found\n')
    end
    
    % Store corners and visualize results (if desired)
    corners(:,:,i) = img_pts;
end

end

