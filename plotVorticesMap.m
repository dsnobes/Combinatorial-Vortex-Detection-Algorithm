function [figs5] = plotVorticesMap(thresholdIntensityVector, velocityData, circleProperties, windingAngleStreamlines, xVortexCore, yVortexCore, radiiVortexCore)
% PLOTVORTICESMAP - Visualizes vortex core locations and vorticity distribution on a map.
%
% Description:
%   This function creates a visual representation of vortices and vorticity data on a map.
%   It highlights the core locations and intensity of the vortices in the given vorticity field.
%
% Inputs:
%   - thresholdIntensityVector: Intensity thresholds for vortex detection.
%   - velocityData: Struct containing velocity field data (rot, x, y).
%   - circleProperties: Properties of circles around vortex cores.
%   - windingAngleStreamlines: Streamlines associated with vortices.
%   - xVortexCore: x-coordinates of vortex cores.
%   - yVortexCore: y-coordinates of vortex cores.
%   - radiiVortexCore: Radii of vortex cores.
%
% Outputs:
%   - figs5: Handle to the generated figure for the vortices map.
%
%--------------------------------------------------------------------------
    
    % Create a new figure
    figs5 = figure(5);
    
    % Apply intensity thresholding to vorticity data
    wthresh = velocityData.rot;
    wthresh(abs(velocityData.rot) < min(thresholdIntensityVector)) = 0;
    wthresh(velocityData.rot > max(thresholdIntensityVector)) = max(thresholdIntensityVector);
    wthresh(velocityData.rot < -max(thresholdIntensityVector)) = -max(thresholdIntensityVector);

    % Plot the vorticity data
    pcolor(velocityData.y, velocityData.x, wthresh);
    shading interp;
    
    % Add a colorbar
    a = colorbar;
    ylabel(a, '\omega (s^{-1})', 'FontSize', 14, 'Rotation', 90, 'FontWeight', 'light', 'FontName', 'Times');
    
    % Set colormap
    jet2 = jet;
    jet2(124:132, :) = ones;
    colormap(jet2); 

    % Hold on to the plot
    hold on;

    % Loop through each vortex
    for vtx = 1:length(windingAngleStreamlines)
        if isempty(circleProperties) == 0
            for k = 1:size(windingAngleStreamlines{vtx}, 2)

                % Plot circles based on the rotation sign of the vortex
                circleProperties(vtx) = circle([xVortexCore(vtx) yVortexCore(vtx)], radiiVortexCore(vtx), 1000, 'k-');
                plot(xVortexCore(vtx), yVortexCore(vtx), 'wx');
            end
        end
    end
    
    % Set axis properties
    axis tight;
end

