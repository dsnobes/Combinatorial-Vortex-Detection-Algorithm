function [drift_vx, drift_vy, driftVelocity, xVortexCore, yVortexCore, radiiVortexCore] = CrossSectionLinesMethod(vortexROIMap, numROIs, velocityData)
% CROSS-SECTIONAL LINES METHOD (CSL)
% Vortex Core Detection and Drift Velocity Calculation

% Description:
%   This function detects vortex cores within regions of interest (ROIs) 
%   representing potential vortices. It uses the Cross-Sectional Lines (CSL)
%   method to pinpoint the core location based on the vorticity distribution.  
%
%   The CSL technique analyzes circumferential velocities along theta and locates
%   the core where the difference between max and min perpendicular velocity 
%   components (v_normal) is greatest. This refined core position enables precise 
%   calculation of vortex drift velocities.
%  
%   The function takes a vortex ROI map and vector field data as input. It outputs
%   the x,y coordinates of vortex cores, core radii, and drift velocities.
%
%   Steps:
%   1) Dilate each ROI to fully capture circumferential velocities
%   2) Extract v_theta distribution within dilated ROI   
%   3) Find max/min v_normal differences along y to locate core row
%   4) Use CSL method to find core x-position within row
%   5) Calculate drift velocities at refined core locations
%
% Inputs:
%   - vortexROIMap: ROI map where each region corresponds to a potential vortex.
%   - numROIs: Total number of detected potential vortices.
%   - velocityData: Updated vector field data.
%
% Outputs:
%   - drift_vx: Drift velocity in the x-direction for each vortex core.
%   - drift_vy: Drift velocity in the y-direction for each vortex core.
%   - driftVelocity: Magnitude of drift velocity for each vortex core.
%   - xVortexCore: x-coordinate of vortex core for each ROI.
%   - yVortexCore: y-coordinate of vortex core for each ROI.
%   - radiiVortexCore: Radius of vortex core for each ROI.
%
% Dependencies:
%   - The function assumes that you have 'strel' and 'imdilate' functions available
%     in your MATLAB environment or as part of external toolboxes.
%
% References:
%   H. Vollmers, Detection of vortices and quantitative evaluation of their 
%   main parameters from experimental velocity data, Meas. Sci. Technol. 12 (2001) 1199. 
%   https://doi.org/10.1088/0957-0233/12/8/329.
%
% Authors: Mathew Bussière, Guilherme Bessa, Bob Koch, and David Nobes.
%          Department of Mechanical Engineering, University of Alberta, 
%          Edmonton, Alberta 
% Contact: dnobes@ualberta.ca
% Version: 1.0
% Date: 10/6/2023
%--------------------------------------------------------------------------

%%

% Initialize variables to store core positions and velocities
core_y_indices = zeros(1, numROIs);
core_x_indices = zeros(1, numROIs);
radiiVortexCore = zeros(1, numROIs);

% Define a structuring element for dilation
dilation_strel = strel('diamond', 1);

% Loop through each region of interest (ROI)
for k = 1:numROIs
    % Dilate the ROI slightly to capture max and min circumferential velocities
    ROI_dilated = vortexROIMap;
    ROI_dilated(vortexROIMap ~= k) = 0;
    
    % Binarize the dilated ROI to create a binary mask
    ROI_binary = ROI_dilated;
    ROI_binary(ROI_binary ~= 0) = 1;
    
    % Dilate the binary ROI to ensure complete velocity capture
    ROI_binary_dilated = imdilate(ROI_binary, dilation_strel);
    
    % Extract velocity components within the ROI
    v_normal = velocityData.vy; 
    v_normal(ROI_binary_dilated ~= 1) = NaN;
    
    % Find indices of max and min velocities along the theta direction
    [maxVelComponent, maxVelComponent_indices] = max(v_normal, [], 2);
    [minVelComponent, minVelComponent_indices] = min(v_normal, [], 2);
    
    % Calculate the difference between max and min velocities to locate the core
    [~, core_y_indices(k)] = max(abs(maxVelComponent - minVelComponent));
    
    % Check for cases where core location is ambiguous or data is missing
    if all(isnan(v_normal(:))) || maxVelComponent(core_y_indices(k)) == minVelComponent(core_y_indices(k)) || abs(minVelComponent_indices(core_y_indices(k)) - maxVelComponent_indices(core_y_indices(k))) < 2
        % Handle special cases where core identification is not possible
        core_y_indices(k) = 1;
        core_x_indices(k) = 1;
        radiiVortexCore(k) = 1;
        minVelComponent(k) = 1;
        maxVelComponent(k) = 1;
    else
        % Calculate x indices corresponding to min and max v_normal values
        core_x_index_min(k) = minVelComponent_indices(core_y_indices(k));
        core_x_index_max(k) = maxVelComponent_indices(core_y_indices(k));
        
        % Swap indices if necessary to ensure correct order
        if core_x_index_min(k) > core_x_index_max(k)
            core_x_index_max_dummy = core_x_index_max(k);
            core_x_index_max(k) = core_x_index_min(k);
            core_x_index_min(k) = core_x_index_max_dummy;
        end
        
        % Calculate vortex core radius based on min and max velocities
        gridSpacing = velocityData.x(end) - velocityData.x(end - 1);
        radiiVortexCore(k) = abs((core_x_index_max(k) - core_x_index_min(k)) * gridSpacing) / 2;
        core_x_center(k) = (maxVelComponent(core_y_indices(k)) + minVelComponent(core_y_indices(k))) / 2;
        
        % Locate x position of the core using CSL method within the specified region
        [~, core_x_indices(k)] = min(abs(v_normal(core_y_indices(k), core_x_index_min(k):core_x_index_max(k)) - core_x_center(k)));
        core_x_indices(k) = core_x_indices(k) + core_x_index_min(k) - 1;
    end
end

% Calculate y and x coordinates of the vortex cores based on indices
yVortexCore = velocityData.y(core_y_indices);
xVortexCore = velocityData.x(core_x_indices); 

% Loop through each core to calculate drift velocities
for k = 1:numel(core_x_indices)
    % Find the closest index in the velocity grid to the core x-coordinate
    [~, x_index(k)] = min(abs(velocityData.x - xVortexCore(k))); 
    % Find the closest index in the velocity grid to the core y-coordinate
    [~, y_index(k)] = min(abs(velocityData.y - yVortexCore(k))); 
    
    % Calculate drift velocity components for each core
    drift_vx(k) = velocityData.vx(y_index(k), x_index(k));
    drift_vy(k) = velocityData.vy(y_index(k), x_index(k));
    driftVelocity(k) = hypot(drift_vx(k), drift_vy(k));
end
