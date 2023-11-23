function figs2 = plotStreamlinesWA(numROIs, xVortexCore, yVortexCore, windingAngleStreamlines, rotationSign, xacc, yacc)
% PLOTSTREAMLINESWA - Visualize streamlines (and their centroids) that meet winding angle criteria.
%
% Description:
%   This function creates a series of figures to visualize streamlines (and their centroids) 
%   that meet winding angle criteria for each vortex core. It also highlights rotation direction 
%   of each vortex.
%
% Inputs:
%   - numROIs: Number of vortex (integer).
%   - xVortexCore: x-coordinate of vortex cores.
%   - yVortexCore: y-coordinate of vortex cores.
%   - windingAngleStreamlines: A cell array of trimmed streamlines that meet winding angle criteria.
%   - rotationSign: Array indicating the rotation sign for each vortex.
%   - xacc: Average of the x-coordinates of the centroids of the streamlines that meet 
%     the wind angle criterion for each vortex.
%   - yacc: Average of the x-coordinates of the centroids of the streamlines that meet 
%     the wind angle criterion for each vortex.
%
% Output:
%   - figs2: Handle to the created figure.
%
%--------------------------------------------------------------------------
    
for vtx = 1:numROIs
 
    figs2 = figure(20 + vtx);
    xlabel('\it y \rm (mm)', 'FontName', 'Times New Roman', 'FontAngle', 'italic');
    ylabel('\it x \rm (mm)', 'FontName', 'Times New Roman', 'FontAngle', 'italic');

    % Set a custom colormap
    colormap([1 1 1]);
    hold on
    
    % Plot streamlines and mark properties
    for k = 1:size(windingAngleStreamlines{vtx}, 2)
        c = windingAngleStreamlines{vtx}{k}; 
        x = c(:, 1);
        y = c(:, 2);
        xa{vtx}(k) = mean(x); % x-coordinates of streamline centroids.
        ya{vtx}(k) = mean(y); % y-coordinates of streamline centroids.

        if rotationSign{vtx} >= 0
            plot(x, y, 'r')
        end
        if rotationSign{vtx} < 0 
            plot(x, y, 'b')
        end
        hold on 
        
        % Plot center of gravity and acceleration point
        plot(xVortexCore(vtx), yVortexCore(vtx), 'k+');
        plot(xa{vtx}, ya{vtx}, 'k.');
        plot(xacc(vtx), yacc(vtx), 'gx');
    end
end
