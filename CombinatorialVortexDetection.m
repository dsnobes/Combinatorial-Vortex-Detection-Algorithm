%--------------------------------------------------------------------------
% Title: COMBINATORIAL VORTEX DETECTION (CVD) ALGORITHM
%
%
% Description:
%   The CVD Algorithm is a research tool designed to streamline and automate vortex analysis in velocity vector fields 
%   obtained from simulations or experiments. Applicable to 2D2C datasets from computational fluid dynamics (CFD) 
%   and particle image velocimetry (PIV), the code implements a processing workflow to identify potential vortices, 
%   confirm and characterize them using advanced methods, and generate visualizations to explore dynamics.
%   The first stage detects possible vortex regions by locating local maxima in vorticity magnitude using 
%   the Maximum Vorticity (MV) method. Around each detected maximum, a Cross-Sectional Lines (CSL) approach 
%   further refines the estimated core location by analyzing circumferential velocities. Regions of Interest (ROIs) 
%   encompassing each potential vortex are then generated.
%   The Winding Angle (WA) method next analyzes the ROIs to definitively identify and characterize vortices based 
%   on a winding angle threshold criterion. Key vortex properties are computed including circulation, area, 
%   core location, rotation direction and drift velocity.
%   Finally, the code produces a serie of visualizations including streamline plots, color maps of vorticity magnitude, 
%   and vorticity profiles. These outputs help researchers visually examine vortex dynamics within the flow, 
%   complementing the quantitative analysis.
%   In summary, the CVD Algorithm automates complex vortex identification and quantification tasks that 
%   previously required manual processing. By streamlining these processes and generating insightful visualizations, 
%   it enables deeper exploration of fluid flows and facilitates research in domains relying on analysis of coherent structures.
%
% Usage:
%   1. Load a vector field ('B00001.vc7') using the 'loadvec' function.
%   2. Adjust the orientation of the vector field for processing. This is a 
%      specific requirement for the current example.
%   3. Set parameters for vortex detection, such as intensity thresholds,
%      structuring elements, and interpolation factor.
%   4. Initialize storage for vortex analysis results.
%   5. Process the vector field frames using the CVD algorithm and store
%      vortex-related data.
%
% Dependencies:
%   - The code assumes that you have the external toolboxes 'Image Processing', 
%     'Signal Processing', 'PIVmat', and 'Wavelet'available in your MATLAB environment.
%
% Authors: Mathew Bussière, Guilherme Bessa, Bob Koch, and David Nobes.
%          Department of Mechanical Engineering, University of Alberta, 
%          Edmonton, Alberta 
% Contact: dnobes@ualberta.ca
% Version: 1.0
% Date: 10/6/2023
%--------------------------------------------------------------------------


%% INPUTS
vectorfield = loadvec('B00001.vc7');
imageInterval = [1]; % Range of image frames to analyze 
timeStep = 2; % Increment between frames to process
timeIncrement = 1.6129e-04; % Time between each frame
        
        %% Adjusting vector field orientation for processing
        %  This is a specific requirement for the current example.
            
            v = vectorfield;
        % Rotate matrices in proper orientation
            v = rotatef(v,pi()/2,'trunc');
        % Move the y-origin and x-origin to the center of the local vector map
            v = shiftf(v,'c');
        % Calculate vorticity and smooth
            w = vec2scal(v, 'curl');
            w = filterf(w, 1, 'gauss', 'same');
        % Transpose matrices for construction of new matrix in proper orientation
            vectorfield.vx = v.vx';
            vectorfield.vy = v.vy';
            vectorfield.rot = w.w';

%% PARAMETERS
thresholdIntensityVector = [0.2, 0.4, 0.6]; % Intensity thresholds for vortex detection
structuringElements = [3, 2, 1]; % Structuring elements for image thresholding

%% INITIALIZE STORAGE 
vortexCirculation = {}; % Circulation for each vortex in each frame
vortexArea = {}; % Area for each vortex in each frame
vortexCircle = {}; % Circle fit parameters for each vortex
vortexLocation = {}; % Vortex core (x,y) coordinates
vortexWAngleStreamlines = {}; % Winding angle for each vortex line 
vortexRotationSign = {}; % Rotation direction (+1/-1) for each vortex
vortexDriftVelocity = {}; % Drift velocity for each vortex
processedFrame = {}; % Track frames processed

%% INICIALIZE COUNTER
counter = 1;

%% MAIN PROCESSING LOOP
tic
for currentImage = imageInterval(1):timeStep:imageInterval(end)
  
  %% Load vector field data
  disp(['Processing frame: ' num2str(currentImage)]);
  vectorFieldData = vectorfield(currentImage);
  fileName = vectorfield(currentImage).name;
  
  %% Detect vortices using Maximum Vorticity method (MV)
  [vortexROIMap, numROIs, velocityData] = MaximumVorticityMethod ...
  (vectorFieldData, thresholdIntensityVector, structuringElements);

  %% Refine vortex locations using Cross-sectional Lines method (CSL)
  [drift_vx, drift_vy, driftVelocity, xVortexCore, yVortexCore, radiiVortexCore] = CrossSectionLinesMethod ...
  (vortexROIMap, numROIs, velocityData);

  %% Generate Regions Of Interest (ROI)
  % Inputs
  skip = 5; % Spacing between streamline points (controls density)
  boxFactorX = 2; % Factor to determine the size of the box around x-direction
  boxFactorY = 2; % Factor to determine the size of the box around y-direction

  [vortexXY, numStreamlines, ROIbox, ROIboxIndex] = RegionsOfInterest ...
  (numROIs, velocityData, xVortexCore, yVortexCore, radiiVortexCore, skip, boxFactorX, boxFactorY, ...
    drift_vx, drift_vy);

  %% Winding Angle method (WA)
  % Inputs
  StartEndDistance = 8; % Threshold distance for start and end points of streamlines
  maxThetaMax = pi()/3; % Maximum angle threshold for winding angle criteria
  thresholdD = 12; % Threshold distance for streamline trimming

  [windingAngleStreamlines, numROIs, bulkStreamlines, xVortexCore, yVortexCore, xacc, yacc,...
      drift_vx, drift_vy, driftVelocity, ROIbox, ROIboxIndex, rotationSign] = WindingAngleMethod ...
  (velocityData, numStreamlines, numROIs, vortexXY, StartEndDistance,...
      maxThetaMax, thresholdD, xVortexCore, yVortexCore, drift_vx, drift_vy,...
      driftVelocity, ROIbox, ROIboxIndex);

  %% Compute circulation and area
  [circulationValues, areaValues, circleProperties, vortexCore] = Circulation ...
  (velocityData, numROIs, bulkStreamlines, xacc, yacc,...
      xVortexCore, yVortexCore, radiiVortexCore);
  
  %% Store data for this frame
  vortexCirculation{counter} = circulationValues;
  vortexArea{counter} = areaValues; 
  vortexCircle{counter} = circleProperties;
  vortexLocation{counter} = vortexCore;
  vortexWAngleStreamlines{counter} = windingAngleStreamlines;
  vortexRotationSign{counter} = rotationSign;
  vortexDriftVelocity{counter} = driftVelocity;
  processedFrame{counter} = currentImage;
  numVortices = numROIs;
  
  %% PLOTS
  
  % Streamlines + Vorticity Field
  [figs1] = plotStreamlines(vortexXY, numStreamlines, ROIboxIndex, numROIs, velocityData, xVortexCore, yVortexCore);

  % Streamlines according to WA:
  [figs2] = plotStreamlinesWA(numROIs, xVortexCore, yVortexCore, windingAngleStreamlines, rotationSign, xacc, yacc);

  % Vorticity and Velocity Profiles
  [figs3, figs4] = plotProfiles(numROIs, velocityData, xVortexCore, yVortexCore, radiiVortexCore, drift_vy, rotationSign);
  
  % Vortices Map
  [figs5] = plotVorticesMap(thresholdIntensityVector, velocityData, circleProperties, windingAngleStreamlines, xVortexCore, yVortexCore, radiiVortexCore);

  %% Increment counter
  counter = counter + 1; 
end
toc

clearvars -except vortexCirculation vortexArea vortexLocation... 
vortexDriftVelocity vortexWAngleStreamlines vortexRotationSign processedFrame numVortices
