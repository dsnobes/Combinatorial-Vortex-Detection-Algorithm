function [vortexROIMap, numROIs, velocityData] = MaximumVorticityMethod(vectorFieldData, ...
    thresholdIntensityVector, structuringElements)
% MAXIMUM VORTICITY METHOD (MV)
%
% Description:
%   This function detects potential vortices in a vector field using the 
%   Maximum Vorticity (MV) method. It is part of the vortex detection workflow
%   in the Combinatorial Vortex Detection (CVD) Algorithm.
%  
%   The MV technique identifies vortex candidates based on local maxima in the 
%   vorticity magnitude. A multilevel thresholding approach is applied using
%   grayscale slicing and morphological operations. This extracts distinct
%   vortical regions at different intensity levels.
%
%   The function takes vector field data and threshold parameters as input. It 
%   outputs a vortex ROI map visualizing candidates, the number of candidates,
%   and the preprocessed vector field data.
%
%   Key steps:
%   1) Apply multilevel thresholding to vorticity magnitude
%   2) Label connected regions at each threshold level 
%   3) Refine regions using morphological opening
%   4) Combine thresholded regions into vortex ROI map
%
% Inputs:
%   - vectorFieldData: A structured dataset containing the vector field components
%     (vx, vy, x, y, rot, mag) to be processed.
%   - thresholdIntensityVector: Intensity thresholds used for vortex detection.
%   - structuringElements: A collection of structuring elements employed for image
%     thresholding operations.
%   - interpolationFactor: Factor used for interpolating velocity vector field data
%     to improve resolution and accuracy.
%
% Outputs:
%   - vortexROIMap: A region-of-interest (ROI) map where each region corresponds to a
%     potential vortex, providing a visual representation of vortex locations.
%   - numROIs: The total number of detected potential vortices, offering a quantitative
%     measure of vortex density.
%   - velocityData: The updated vector field data after preprocessing and analysis.
%
% References:
%   R.C. Strawn, D.N. Kenwright, J. Ahmad, Computer Visualization of Vortex 
%   Wake Systems, AIAA Journal. 37 (1999) 511–512. https://doi.org/10.2514/2.744.
%
% Authors: Mathew Bussière, Guilherme Bessa, Bob Koch, and David Nobes.
%          Department of Mechanical Engineering, University of Alberta, 
%          Edmonton, Alberta 
% Contact: dnobes@ualberta.ca
% Version: 1.0
% Date: 10/6/2023
%--------------------------------------------------------------------------

%%

% Rename the input for clarity
inputData = vectorFieldData;

% find the indices of all non-zero elements
[row_indices, col_indices] = find(inputData.vx);
x0  = min(row_indices); xend = max(row_indices);
y0  = min(col_indices); yend = max(col_indices);

% Exctracting the non-zero elements
inputData.x = inputData.x(x0:xend);
inputData.y = inputData.y(y0:yend);
inputData.vx = inputData.vx(x0:xend,y0:yend);
inputData.vy = inputData.vy(x0:xend,y0:yend);
inputData.mag = hypot(inputData.vx, inputData.vy);
inputData.rot = inputData.rot(x0:xend,y0:yend);

% Multilevel threshold
% Define structuring elements
s1 = strel('diamond', structuringElements(1));
s2 = strel('diamond', structuringElements(2));
s3 = strel('diamond', structuringElements(3));

% Define the threshold intensity vector (vv)
vv = thresholdIntensityVector;

%% POSITIVE THRESHOLD

vvp = vv;  % Assign the threshold intensity vector to vvp

% Apply grayscale slicing to the rotation data using the threshold vector
Xp = grayslice(inputData.rot, vvp);
Xp = double(Xp);  % Convert the thresholded data to double precision for computations
Xp_binary = Xp;  % Create a copy of the thresholded data

% Convert non-zero elements to 1, creating a binary mask
Xp_binary(Xp ~= 0) = 1;

% Apply morphological opening using structuring element s1
Xp_binary = imopen(Xp_binary, s1);

% Label connected components in the binary image
[L, num] = bwlabel(Xp_binary, 4);

% Preallocate variables for positive thresholding
xp1 = cell(1, num);
xp1t = zeros(size(L));
xp11 = cell(1, num);
xp11_binary = cell(1, num);
xp2 = cell(1, num);
xp2t = zeros(size(L));
xp22 = cell(1, num);
xp22_binary = cell(1, num);
xp3 = cell(1, num);
xp3t = zeros(size(L));

% Preallocate variables for labeling
L1 = cell(1, num);
num1 = zeros(1, num);
L2 = cell(1, num);
num2 = cell(1, num);

% Process each labeled region from the first threshold
for i = 1:num
    % Create a copy of the thresholded data for the current label
    xp1{i} = Xp;
    % Zero out regions not belonging to the current label
    xp1{i}(L ~= i) = 0;

    % Peel back the first layer to expose the groups inside each label of the first threshold
    xp11{i} = xp1{i};
    xp11{i}(xp1{i} == 1) = 0;

    % Apply morphological opening to the peeled layer using structuring element s2
    xp11_binary{i} = xp11{i};
    xp11_binary{i}(xp11{i} ~= 0) = 1;
    xp11_binary{i} = imopen(xp11_binary{i}, s2);

    % Label connected components in the peeled layer
    [L1{i}, num1(i)] = bwlabel(xp11_binary{i}, 4);
    
    if num1(i) > 0
        % Process each labeled region from the second threshold
        for j = 1:num1(i)
            % Create a copy of the peeled data for the current label
            xp2{i}{j} = xp11{i};
            % Zero out regions not belonging to the current label
            xp2{i}{j}(L1{i} ~= j) = 0;

            % Peel back the second layer to expose the groups inside each label of the second threshold
            xp22{i}{j} = xp2{i}{j};
            xp22{i}{j}(xp2{i}{j} == 2) = 0;

            % Apply morphological opening to the second peeled layer using structuring element s3
            xp22_binary{i}{j} = xp22{i}{j};
            xp22_binary{i}{j}(xp22{i}{j} ~= 0) = 1;
            xp22_binary{i}{j} = imopen(xp22_binary{i}{j}, s3);

            % Label connected components in the second peeled layer
            [L2{i}{j}, num2{i}{j}] = bwlabel(xp22_binary{i}{j}, 4);

            if num2{i}{j} >= 2
                % Process each labeled region from the third threshold
                for k = 1:num2{i}{j}
                    % Create a copy of the second peeled data for the current label
                    xp3{i}{j}{k} = xp22{i}{j};
                    % Zero out regions not belonging to the current label
                    xp3{i}{j}{k}(L2{i}{j} ~= k) = 0;

                    % Clear intermediate layers to avoid overlapping regions
                    xp1{i} = zeros(size(xp1{i}));
                    xp2{i}{j} = zeros(size(xp2{i}{j}));

                    % Add the third thresholded region to the cumulative result
                    xp3{i}{j}{k}(xp3{i}{j}{k} ~= 0) = 1;
                    xp3t = xp3t + xp3{i}{j}{k};
                end
            else
                % Clear intermediate layers
                if num1(i) >= 2
                    xp1{i} = zeros(size(xp1{i}));
                end
                if num1(i) < 2
                    xp2{i}{j} = zeros(size(xp2{i}{j}));
                end
            end

            % Add the second thresholded region to the cumulative result
            xp2{i}{j}(xp2{i}{j} ~= 0) = 1;
            xp2t = xp2t + xp2{i}{j};
        end
    end

    % Add the first thresholded region to the cumulative result
    xp1{i}(xp1{i} ~= 0) = 1;
    xp1t = xp1t + xp1{i};
end

% Combine thresholded regions from all levels
xptt = xp1t + xp2t + xp3t;

% Label positive ROIs
[ROI_pos, roipos_num] = bwlabel(xptt, 4);

% Clear intermediate variables
clear L L1 L2 num num1 num2

%% NEGATIVE THRESHOLD

% Reverse the threshold intensity vector for negative thresholding
vvn = wrev(-vv);

% Apply grayscale slicing to the rotation data using the reversed threshold vector
Xn = grayslice(inputData.rot, vvn);

% Convert the thresholded data to double precision for computations
Xn = double(Xn);

% Invert the data and add a constant to enhance visibility
Xn = -Xn + 3;

% Create a copy of the inverted data for further processing
Xn_binary = Xn;

% Convert non-zero elements to 1, creating a binary mask
Xn_binary(Xn ~= 0) = 1;

% Apply morphological opening using structuring element s1
Xn_binary = imopen(Xn_binary, s1);

% Label connected components in the binary image
[L, num] = bwlabel(Xn_binary, 4);

% Preallocate variables for negative thresholding
xn1 = cell(1, num);
xn1t = zeros(size(L));
xn11 = cell(1, num);
xn11_binary = cell(1, num);
xn2 = cell(1, num);
xn2t = zeros(size(L));
xn22 = cell(1, num);
xn22_binary = cell(1, num);
xn3 = cell(1, num);
xn3t = zeros(size(L));

% Preallocate variables for labeling
L1 = cell(1, num);
num1 = zeros(1, num);
L2 = cell(1, num);
num2 = cell(1, num);

% Process each labeled region from the first threshold
for i = 1:num
    % Create a copy of the thresholded data for the current label
    xn1{i} = Xn;
    % Zero out regions not belonging to the current label
    xn1{i}(L ~= i) = 0;

    % Peel back the first layer to expose the groups inside each label of the first threshold
    xn11{i} = xn1{i};
    xn11{i}(xn1{i} == 1) = 0;

    % Apply morphological opening to the peeled layer using structuring element s2
    xn11_binary{i} = xn11{i};
    xn11_binary{i}(xn11{i} ~= 0) = 1;
    xn11_binary{i} = imopen(xn11_binary{i}, s2);

    % Label connected components in the peeled layer
    [L1{i}, num1(i)] = bwlabel(xn11_binary{i}, 4);

    if num1(i) > 0
        % Process each labeled region from the second threshold
        for j = 1:num1(i)
            % Create a copy of the peeled data for the current label
            xn2{i}{j} = xn11{i};
            % Zero out regions not belonging to the current label
            xn2{i}{j}(L1{i} ~= j) = 0;

            % Peel back the second layer to expose the groups inside each label of the second threshold
            xn22{i}{j} = xn2{i}{j};
            xn22{i}{j}(xn2{i}{j} == 2) = 0;

            % Apply morphological opening to the second peeled layer using structuring element s3
            xn22_binary{i}{j} = xn22{i}{j};
            xn22_binary{i}{j}(xn22{i}{j} ~= 0) = 1;
            xn22_binary{i}{j} = imopen(xn22_binary{i}{j}, s3);

            % Label connected components in the second peeled layer
            [L2{i}{j}, num2{i}{j}] = bwlabel(xn22_binary{i}{j}, 4);

            if num2{i}{j} >= 2
                % Process each labeled region from the third threshold
                for k = 1:num2{i}{j}
                    % Create a copy of the second peeled data for the current label
                    xn3{i}{j}{k} = xn22{i}{j};
                    % Zero out regions not belonging to the current label
                    xn3{i}{j}{k}(L2{i}{j} ~= k) = 0;

                    % Clear intermediate layers to avoid overlapping regions
                    xn1{i} = zeros(size(xn1{i}));
                    xn2{i}{j} = zeros(size(xn2{i}{j}));

                    % Add the third thresholded region to the cumulative result
                    xn3{i}{j}{k}(xn3{i}{j}{k} ~= 0) = 1;
                    xn3t = xn3t + xn3{i}{j}{k};
                end
            else
                % Clear intermediate layers
                if num1(i) >= 2
                    xn1{i} = zeros(size(xn1{i}));
                end
                if num1(i) < 2
                    xn2{i}{j} = zeros(size(xn2{i}{j}));
                end
            end

            % Add the second thresholded region to the cumulative result
            xn2{i}{j}(xn2{i}{j} ~= 0) = 1;
            xn2t = xn2t + xn2{i}{j};
        end
    end

    % Add the first thresholded region to the cumulative result
    xn1{i}(xn1{i} ~= 0) = 1;
    xn1t = xn1t + xn1{i};
end

% Combine thresholded regions from all levels
xntt = xn1t + xn2t + xn3t;

% Label negative ROIs
[ROI_neg, roineg_num] = bwlabel(xntt, 4);


%% Combine the positive and negative ROIs into one ROI map

% Increment the labels of negative ROIs by the number of positive ROIs
% This ensures that the labels of the combined ROI map don't overlap between positive and negative ROIs
ROI_neg(ROI_neg ~= 0) = ROI_neg(ROI_neg ~= 0) + roipos_num;

% Set the labels of positive ROIs in the combined ROI map to 0
% This removes the positive ROIs from the combined map since their labels are no longer needed
ROI_pos(ROI_neg ~= 0) = 0;

% Sum the positive and negative ROIs to create the final combined ROI map
vortexROIMap = ROI_neg + ROI_pos;

% Calculate the total number of ROIs by summing the numbers of positive and negative ROIs
numROIs = roipos_num + roineg_num;

% Convert the combined ROI map to double precision for consistency
vortexROIMap = double(vortexROIMap);

%% Display result
figure(4)
pcolor(vortexROIMap)

%% Rename the output for clarity
velocityData = inputData;

%% Clear unnecessary variables
clearvars -except vortexROIMap numROIs velocityData

end
