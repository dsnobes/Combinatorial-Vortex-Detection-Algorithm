function [vortexXY, numStreamlines, ROIbox, ROIboxIndex] = RegionsOfInterest(numROIs, velocityData, xVortexCore, yVortexCore, radiiVortexCore, skip, boxFactorX, boxFactorY, drift_vx, drift_vy)
% REGIONS OF INTEREST (ROI)
% Function for identifying vortex core bounding boxes and extracting streamlines within them.

% Description:
%   This function defines regions of interest (ROIs) around potential vortex cores 
%   detected in prior steps. It extracts streamlines seeded within each ROI 
%   to analyze the local flow behavior.  
%
%   Bounding boxes around each vortex core are defined based on the core's 
%   location, radius, and velocity. The box dimensions aim to fully encompass  
%   the vortex region.
%
%   Velocity data is extracted within each ROI box. Streamlines are then
%   generated with regular starting point spacing to characterize the flow.
%
%   The number of streamlines in each vortex provides a metric of complexity.
%   The streamline trajectories and ROI boxes enable detailed vortex visualization.
%

% Inputs:
%   - numROIs: The number of vortex cores (integer).
%   - velocityData: A structured dataset containing velocity field data.
%   - xVortexCore: The x-coordinate of vortex cores.
%   - yVortexCore: The y-coordinate of vortex cores.
%   - radiiVortexCore: The radius of vortex cores.
%   - skip: The spacing between streamline points (controls streamline density).
%   - boxFactorX: A factor used to determine the size of the box around the x-direction of each core.
%   - boxFactorY: A factor used to determine the size of the box around the y-direction of each core.
%   - drift_vx: Drift velocity in the x-direction for each vortex core.
%   - drift_vy: Drift velocity in the y-direction for each vortex core.

% Outputs:
%   - vortexXY: A cell array containing streamlines for each vortex core, providing a
%     detailed representation of the fluid flow behavior within the cores.
%   - numStreamlines: An array indicating the number of streamlines for each vortex core,
%     which is a key metric for assessing vortex core complexity.
%   - ROIbox: A cell array describing the bounding boxes around vortex cores,
%     allowing the visualization of the core regions of interest.
%   - ROIboxIndex: A cell array indicating the indices of bounding box coordinates for
%     further data extraction and analysis.
%
% Dependencies:
%   - This function doesn't depend on external MATLAB toolboxes or functions.
%
%
% Authors: Mathew Bussière, Guilherme Bessa, Bob Koch, and David Nobes.
%          Department of Mechanical Engineering, University of Alberta, 
%          Edmonton, Alberta 
% Contact: dnobes@ualberta.ca
% Version: 1.0
% Date: 10/6/2023
%--------------------------------------------------------------------------


%%

% Initialize variables
clear n m X Y Xs Ys sx sy
numStreamlines = zeros(1, numROIs);
vortexXY = cell(1, numROIs);
vvxBox = cell(1, numROIs);
vvyBox = cell(1, numROIs);
vxBox = cell(1, numROIs);
vyBox = cell(1, numROIs);
ROIbox = cell(1, numROIs);
ROIboxIndex = cell(1, numROIs);

% Loop through each vortex core
for vtx = 1:numROIs
    ref_velx{vtx} = velocityData.vx - drift_vx(vtx);
    ref_vely{vtx} = velocityData.vy - drift_vy(vtx);
    
    % Define bounding box around the vortex core
    ROIbox{vtx} = [xVortexCore(vtx) - boxFactorX * radiiVortexCore(vtx), xVortexCore(vtx) + boxFactorX * radiiVortexCore(vtx), yVortexCore(vtx) - boxFactorY * radiiVortexCore(vtx), yVortexCore(vtx) + boxFactorY * radiiVortexCore(vtx)];
    
    % Calculate indices for the bounding box coordinates
    [~, xmin_index(vtx)] = min(abs(velocityData.x - ROIbox{vtx}(1)));
    [~, xmax_index(vtx)] = min(abs(velocityData.x - ROIbox{vtx}(2)));
    [~, ymin_index(vtx)] = min(abs(velocityData.y - ROIbox{vtx}(3)));
    [~, ymax_index(vtx)] = min(abs(velocityData.y - ROIbox{vtx}(4)));
    ROIboxIndex{vtx} = [xmin_index(vtx), xmax_index(vtx), ymin_index(vtx), ymax_index(vtx)];

    % Extract velocity components within the bounding box
    vvxBox{vtx} = ref_velx{vtx}(ROIboxIndex{vtx}(3):ROIboxIndex{vtx}(4), ROIboxIndex{vtx}(1):ROIboxIndex{vtx}(2));
    vvyBox{vtx} = ref_vely{vtx}(ROIboxIndex{vtx}(3):ROIboxIndex{vtx}(4), ROIboxIndex{vtx}(1):ROIboxIndex{vtx}(2));
    vxBox{vtx} = velocityData.x(ROIboxIndex{vtx}(1):ROIboxIndex{vtx}(2));
    vyBox{vtx} = velocityData.y(ROIboxIndex{vtx}(3):ROIboxIndex{vtx}(4));

    % Specify grid points for streamline starting points
    [n, m] = size(vvxBox{vtx});
    [X, Y] = meshgrid(vxBox{vtx}, vyBox{vtx});
    [Xs, Ys] = meshgrid(vxBox{vtx}(1):skip:vxBox{vtx}(end), vyBox{vtx}(1):skip:vyBox{vtx}(end));

    % Reshape grid points into vectors for starting points
    sx = reshape(Xs, 1, numel(Xs));
    sy = reshape(Ys, 1, numel(Ys));
    numStreamlines(vtx) = numel(sy);

    % Generate streamlines within the vortex core bounding box
    clear XYd
    XYd = stream2(X, Y, vvxBox{vtx}, vvyBox{vtx}, sx, sy, [0.25, 86]);
    vortexXY{vtx} = XYd;
end
