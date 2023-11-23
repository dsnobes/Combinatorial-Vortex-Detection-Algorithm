function [circulationValues, areaValues, circleProperties, vortexCore] =...
 Circulation(velocityData, numROIs, bulkStreamlines, xacc, yacc, xVortexCore, yVortexCore, radiiVortexCore)
% Calculate Core-Region Circulation of Detected Vortices
%
% Description:
%   This function computes the circulation and area of the core region for
%   detected vortices. The function processes each detected vortex and calculates the circulation and area of
%   its core region using the provided streamline data. It fits a circle to the core region,
%   identifies grid points inside the region defined by the circle, and computes the circulation
%   and area based on vorticity values within the core region. The core region's centroid coordinates
%   and fitted circle properties are also determined. If no streamline data is available for a vortex, 
%   the corresponding outputs for that vortex are empty.
%
% Inputs:
%   - velocityData: Structure containing velocity field data.
%   - numROIs: Number of detected vortex cores (integer).
%   - bulkStreamlines: Cell array containing concatenated streamline points for each vortex core.
%   - xacc: Array of x-coordinates of streamline centroids.
%   - yacc: Array of y-coordinates of streamline centroids.
%   - xVortexCore: Array of x-coordinates of vortex cores.
%   - yVortexCore: Array of y-coordinates of vortex cores.
%   - radiiVortexCore: Array of radii for vortex cores.
%
% Outputs:
%   - circulationValues: Array of core-region circulation values for each vortex core.
%   - areaValues: Array of core-region area values for each vortex core.
%   - circleProperties: Array of properties defining the fitted core region circles.
%   - vortexCore: Array of coordinates of the vortex core centroids.
% Dependencies:
% - inpolygon: Function to determine if points are inside or on edge of a polygonal region (MATLAB-defined)
% - circle: Function to fit a circle to streamline points (custom-defined). 
%
% Authors: Mathew Bussière, Guilherme Bessa, Bob Koch, and David Nobes.
%          Department of Mechanical Engineering, University of Alberta, 
%          Edmonton, Alberta 
% Contact: dnobes@ualberta.ca
% Version: 1.0
% Date: 10/6/2023
%%

% Set up the grid coordinates
x = velocityData.y; 
y = velocityData.x;

% Calculate the area of a single pixel assuming equal dimensions
dA = abs((x(2)-x(1))*(y(2)-y(1)));

% Create a meshgrid for the entire domain
[X,Y] = meshgrid(x, y);
x = reshape(X, 1, numel(X))';
y = reshape(Y, 1, numel(Y))';

% Initialize output variables
circulationValues = zeros(1, length(bulkStreamlines));
areaValues = zeros(1, length(bulkStreamlines));
circleProperties = zeros(1, length(bulkStreamlines));
vortexCore = zeros(2, length(bulkStreamlines));
vortexLocations = cell(1, length(bulkStreamlines));

if numROIs > 1
    % Loop through each vortex
    for vtx = 1:numROIs
        if isempty(bulkStreamlines{vtx}) == 1
            % Skip if no streamline data for this vortex
            circulationValues(vtx) =[]; 
            areaValues(vtx) = [];
            circleProperties(vtx) = []; 
            vortexCore(:,vtx) = [];
            vortexLocations{vtx} = [];
        else
            % Calculate circulation and area for vortex with streamline data

            % Fit an circle to the streamline points to define the core region
            circleProperties(vtx) = circle([xVortexCore(vtx),yVortexCore(vtx)],radiiVortexCore(vtx),1000,'-');
            XV = get(circleProperties(vtx),'XData');
            YV = get(circleProperties(vtx),'YData');

            % Use inpolygon to identify grid points inside the ellipse
            in = inpolygon(X, Y, XV, YV);
            in = double(in);
            IN{vtx} = in;
            C = in;

            % Store only the vorticity that lies within the ellipse in variable C 
            C(in == 1) = velocityData.rot(in == 1);

            % Remove zeros and calculate circulation and area
            C(C == 0) = [];
            circulationValues(vtx) = dA * sum(C);
            areaValues(vtx) = dA * length(C);

            % Store vortex core coordinates
            vortexCore(1, vtx) = xacc(vtx);
            vortexCore(2, vtx) = yacc(vtx);

            vortexLocations{vtx} = [vortexCore(1, vtx), vortexCore(2, vtx)];
        end      
    end
else
    % If only one vortex or no vortex
    circulationValues =[]; 
    areaValues = [];
    circleProperties = []; 
    vortexCore= [];
    vortexLocations= [];
end
