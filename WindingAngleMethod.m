function [windingAngleStreamlines, numROIs, bulkStreamlines, xVortexCore, yVortexCore, xacc, yacc, drift_vx, drift_vy, driftVelocity, ROIbox, ROIboxIndex, rotationSign] =...
 WindingAngleMethod(velocityData, numStreamlines, numROIs, vortexXY, StartEndDistance, maxThetaMax, thresholdD, xVortexCore, yVortexCore, drift_vx, drift_vy, driftVelocity, ROIbox, ROIboxIndex)
% WINDING ANGLE METHOD (WA)
% Analyzing Winding Angles and Streamline Criteria for Vortex Core Detection

% Description:
%   This function performs vortex identification and characterization on a 
%   velocity field using the winding angle method. It is part of the vortex  
%   detection workflow in the Combinatorial Vortex Detection (CVD) Algorithm.
% 
%   The function takes as input velocity data, number of streamlines, vortex 
%   candidate locations from prior steps, and parameters for winding angle  
%   thresholds. It traces streamlines originating from the vortex candidates,
%   computes the winding angle for each streamline, and retains candidates
%   meeting winding angle threshold criteria.
%
%   For confirmed vortices, it calculates key properties including circulation,
%   area, core location, and drift velocity. It also refines the vortex core
%   position using the streamline centroids. Useful visualizations like 
%   streamtraces, vorticity maps, and vortex boundaries are generated.
%
%   Outputs include the refined vortex core locations, magnitudes, signs of 
%   rotation, streamlines, and vortex region boundaries. By 
%   automating the winding angle analysis, this function enables robust vortex
%   identification and characterization.

% Inputs:
% - velocityData: A structure containing velocity field data.
% - numStreamlines: An array indicating the number of streamlines for each vortex candidate.
% - numROIs: The total number of candidates vortices (integer).
% - vortexXY: A cell array containing streamlines for each vortex candidate.
% - StartEndDistance: A threshold distance for start and end points of streamlines.
% - maxThetaMax: The maximum angle threshold for winding angle criteria.
% - thresholdD: A threshold distance for streamline trimming.
% - xVortexCore: The x-coordinate of the core of each vortex candidate.
% - yVortexCore: The y-coordinate of the core of each vortex candidate..
% - drift_vx: An array of x-component of drift velocity.
% - drift_vy: An array of y-component of drift velocity.
% - driftVelocity: An array of drift velocity magnitudes.
% - ROIbox: A cell array describing bounding boxes around vortex candidate.
% - ROIboxIndex: A cell array indicating indices of bounding box coordinates.

% Outputs:
% - windingAngleStreamlines: A cell array of trimmed streamlines that meet winding angle criteria.
% - numROIs: The updated number of vortices after analysis.
% - bulkStreamlines: A cell array of concatenated streamline points for vortex cores.
% - xVortexCore: The updated x-coordinates of vortex cores.
% - yVortexCore: The updated y-coordinates of vortex cores.
% - xacc: An array of x-coordinates of streamline centroids.
% - yacc: An array of y-coordinates of streamline centroids.
% - drift_vx: The updated array of x-component of drift velocity.
% - drift_vy: The updated array of y-component of drift velocity.
% - driftVelocity: The updated array of drift velocity magnitudes.
% - ROIbox: The updated cell array of bounding boxes around vortex cores.
% - ROIboxIndex: The updated cell array of indices of bounding box coordinates.
% - rotationSign: A cell array of signs indicating the direction of rotation of the vortices.
%
% Dependencies:
%   - This function doesn't depend on external MATLAB toolboxes or functions.
%
% References:
%   L.M. Portela, Identification and characterization of vortices in the 
%   turbulent boundary layer, Ph.D., Stanford Uni-versity, 1998.
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
bulkStreamlines = {};
vortexcontour = {};

% Calculate alpha
for vtx = 1:numROIs
    % Clear variables
    clear theta_wk
    clear theta
    clear angle 
    clear SLw
    clear thetasum
    clear thetasum_tot
    clear xw
    clear yw
    clear c
    clear k
    clear xynorm
    clear n
    
    % Preallocate variables
    SL_size = size(vortexXY{vtx}, 2);
    theta = zeros(1, SL_size);
    yw = zeros(1, SL_size);
    xw = zeros(1, SL_size);
    
    % Calculate the angles between grid-points for all streamlines
    for k = 1:numStreamlines(vtx)
        c = vortexXY{vtx}{k};
        if isempty(c) == 0
            streamlines{vtx}(k) = k;
            for i = 1:size(c, 1)
                xw(i) = c(i, 1);
                yw(i) = c(i, 2);

                theta(k, i) = 0;
                trig = 0;
                if i > 2
                    v1 = [(xw(i - 1) - xw(i - 2)), (yw(i - 1) - yw(i - 2)), 0];
                    v2 = [(xw(i) - xw(i - 1)), (yw(i) - yw(i - 1)), 0];
                    vn = [0, 0, 1];

                    theta(k, i) = acos((dot(v1, v2)) / ((norm(v1) * norm(v2))));
                    sign = dot(vn, cross(v1, v2));

                    if sign < 0
                        theta(k, i) = -theta(k, i);
                    end 
                    thetasum{k}(i) = sum(theta(k, :));
                    trig = k;
                end

            end

            half(k) = round(numel(theta(k, :)) / 2);
            fh(k) = abs(sum(theta(k, 1:half(k))));
            lh(k) = abs(sum(theta(k, half(k):end)));
            
            if fh(k) < lh(k)
                theta(k, :) = fliplr(theta(k, :));
                for i = 1:size(c, 1)
                    thetasum{k}(i) = sum(theta(k, 1:i));
                end
                inv(k) = 1;
            else
                inv(k) = 0;
            end
            
            if trig == k
                if thetasum{k}(end) >= 0
                    n(k) = (floor(abs(thetasum{k}(end) / (1 * pi()))));
                else
                    n(k) = -(floor(abs(thetasum{k}(end) / (1 * pi()))));
                end   
                cutoff(k) = n(k) * 1 * pi();
                [min_difference(k), cut_index(k)] = min(abs(thetasum{k}(:) - cutoff(k)));
            else
                n(k) = 0;
                cut_index(k) = 1;
                thetasum{k} = 0;
            end
            maxtheta = max(abs(theta), [], 2);
        end
    end

n = n';
n(isnan(n) == 1) = 0;

% Rebuild all streamlines with the chopped segments deleted
for k = 1:numStreamlines(vtx)
    if inv(k) == 1
        d = flipud(vortexXY{vtx}{k});
    else
        d = (vortexXY{vtx}{k});
    end
    
    xwa = d(1:cut_index(k), 1);
    ywa = d(1:cut_index(k), 2);
    xywa = [xwa ywa];
    SLc{k} = xywa;
    
    xlength = xwa(end) - xwa(1);
    ylength = ywa(end) - ywa(1);
    xynorm(k) = sqrt(xlength^2 + ylength^2);
    
    clear xwa
    clear ywa
    clear xywa
    clear d
end
xynorm = xynorm';

% Total sum of the angles at the end of the streamline
clear theta_wk
theta_wk = n;
Theta{vtx} = n;
SLw = streamlines{vtx}(:);
theta_wk((abs(n) == 0) | (xynorm > StartEndDistance) | (maxtheta >= maxThetaMax)) = [];
% rotationSign{vtx} = theta_wk;

% Identify streamlines that fit winding angle criteria and label them as SLw
SLw((abs(n) == 0) | (xynorm > StartEndDistance) | (maxtheta >= maxThetaMax)) = [];
if isempty(SLw) == 0
    for k = 1:size(SLw, 1)
        d = (SLc{SLw(k)});
        xwa = d(:, 1);
        ywa = d(:, 2);
        xywa = [xwa ywa];
        windingAngleStreamlines{vtx}{k} = xywa;
        
        clear xwa
        clear ywa
        clear xywa
        clear d
    end
    
    for k = 1:size(SLw, 1)
        c = windingAngleStreamlines{vtx}{k};
        x = c(:, 1);
        y = c(:, 2);
        xa{vtx}(k) = mean(x);
        ya{vtx}(k) = mean(y);
        
        if theta_wk(k) >= 0
            rotationSign{vtx} = 1;
        end
        
        if theta_wk(k) < 0 
            rotationSign{vtx} = -1;
        end
    end
    
    vnum{vtx} = zeros(1, size(xa{vtx}, 2));
    vnum{vtx}(1) = 1;
    
    if size(xa{vtx}, 2) > 1
        for k = 2:size(xa{vtx}, 2)
            dc{vtx}(k) = sqrt(((xa{vtx}(k) - xa{vtx}(1))^2) + ((ya{vtx}(k) - ya{vtx}(1))^2));
        end

        windingAngleStreamlines{vtx}(dc{vtx} > thresholdD) = [];
        vnum{vtx}(dc{vtx} > thresholdD) = [];
        xa{vtx}(dc{vtx} > thresholdD) = [];
        ya{vtx}(dc{vtx} > thresholdD) = [];
    else
        windingAngleStreamlines{vtx} = [];
        vnum{vtx} = [];
        xa{vtx} = [];
        ya{vtx} = [];
    end 
else
    windingAngleStreamlines{vtx} = [];
    vnum{vtx} = [];
    rotationSign{vtx} = [];
end

disp(vtx)
end

for vtx = 1:length(windingAngleStreamlines)
    if isempty(windingAngleStreamlines{vtx}) == 0
        xacc(vtx) = mean(xa{vtx});
        yacc(vtx) = mean(ya{vtx});
    else
        xacc(vtx) = NaN;
        yacc(vtx) = NaN;
    end
end

% For boxes with more than one set of closed streamlines, keep the streamlines
% belonging to the largest set
clear SLwa_keep
for vtx = 1:numROIs
    SLwa_keep{vtx} = find(vnum{vtx} == mode(vnum{vtx}));
    windingAngleStreamlines{vtx} = windingAngleStreamlines{vtx}(SLwa_keep{vtx});
    vnum{vtx} = vnum{vtx}(SLwa_keep{vtx});
    
    X = [];
    Y = [];
    if isempty(windingAngleStreamlines{vtx}) == 0
        for k = 1:length(windingAngleStreamlines{vtx})
            cnt = windingAngleStreamlines{vtx}{k};
            X = [X; cnt(:, 1)];
            Y = [Y; cnt(:, 2)];
            bulkStreamlines{vtx} = [X Y];
            vortexcontour{vtx} = boundary(bulkStreamlines{vtx}(:, 1), bulkStreamlines{vtx}(:, 2), 0.2);
        end
    else
        bulkStreamlines{vtx} = [];
        vortexcontour{vtx} = [];
    end
end
clear cnt

% Remove ROIs that do not satisfy the winding angle criteria
notvtx = [];
i = 0;
for vtx = 1:numROIs
    if isempty(windingAngleStreamlines{vtx})
        i = i + 1;
        xVortexCore(vtx) = NaN;
        yVortexCore(vtx) = NaN;
        drift_vy(vtx) = NaN;
        drift_vx(vtx) = NaN;
        dirr(i) = vtx;
        notvtx = [notvtx vtx];
    else
        [~, idx_xcore] = find(velocityData.x == xVortexCore(vtx));
        [~, idx_ycore] = find(velocityData.y == yVortexCore(vtx));
        [~, idx_xacc] = min(abs(velocityData.y - xacc(vtx)));
        [~, idx_yacc] = min(abs(velocityData.x - yacc(vtx)));
    
        if velocityData.rot(idx_ycore, idx_xcore) < 0 && velocityData.rot(idx_yacc, idx_xacc) >= 0
            notvtx = [notvtx vtx];
        elseif velocityData.rot(idx_ycore, idx_xcore) >= 0 && velocityData.rot(idx_yacc, idx_xacc) < 0
            notvtx = [notvtx vtx];
        end
    end
    clear idx_xcore idx_ycore idx_xacc idx_yacc
end

if i ~= 0
    rotationSign(dirr) = [];
end

xacc(notvtx) = [];
yacc(notvtx) = [];
xVortexCore(notvtx) = [];
yVortexCore(notvtx) = [];
drift_vx(notvtx) = [];
drift_vy(notvtx) = [];
driftVelocity(notvtx) = [];
streamlines(notvtx) = [];
windingAngleStreamlines(notvtx) = [];
numROIs = length(windingAngleStreamlines);

ROIboxIndex(notvtx) = [];
vortexcontour(notvtx) = [];
bulkStreamlines(notvtx) = [];
ROIbox(notvtx) = [];
