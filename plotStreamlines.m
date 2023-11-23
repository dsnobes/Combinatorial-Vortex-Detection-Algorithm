function [figs1] = plotStreamlines(vortexXY, numStreamlines, ROIboxIndex, numROIs, velocityData, xVortexCore, yVortexCore)
% PLOTSTREAMLINES - Visualize streamlines within vortex regions and associated vorticity magnitude.
%
% Description:
%   This function creates a series of figures to visualize streamlines within vortex regions
%   alongside the vorticity magnitude as a background image.
%
% Inputs:
%   - vortexXY: Cell array containing streamlines for each vortex core.
%   - numStreamlines: Array indicating the number of streamlines for each vortex core.
%   - ROIboxIndex: Cell array indicating indices of bounding box coordinates.
%   - numROIs: Number of vortex cores (integer).
%   - velocityData: Structure containing velocity field data.
%   - xVortexCore: x-coordinate of vortex cores.
%   - yVortexCore: y-coordinate of vortex cores.
%
% Output:
%   - figs1: Handle to the created figure.
%
%--------------------------------------------------------------------------


for vtx = 1:numROIs
    % Create a figure with the velocity magnitude as the background image
    % along with the computed streamlines
    figs1 = figure(10 + vtx);
    x_min = velocityData.y(ROIboxIndex{vtx}(1));
    x_max = velocityData.y(ROIboxIndex{vtx}(2));
    y_min = velocityData.x(ROIboxIndex{vtx}(3));
    y_max = velocityData.x(ROIboxIndex{vtx}(4));
    rot_box{vtx} = velocityData.rot(ROIboxIndex{vtx}(3):ROIboxIndex{vtx}(4), ROIboxIndex{vtx}(1):ROIboxIndex{vtx}(2));
    
    % Display velocity magnitude as the background image
    imagesc([x_min x_max], [y_min y_max], rot_box{vtx});
    xlabel('\it x \rm (mm)', 'FontName', 'Times New Roman', 'FontAngle', 'italic');
    ylabel('\it y \rm (mm)', 'FontName', 'Times New Roman', 'FontAngle', 'italic');
    % Set colormap
        %jet2 = jet;
        %jet2(124:132,:) = ones;
    colormap(jet)
    caxis([-0.5 1])
    hold on
    
    % Plot streamlines
    for k = 1:1:numStreamlines(vtx)
        c = vortexXY{vtx}{k};
        if isempty(c) == 0
            x = c(:, 1);
            y = c(:, 2);
            plot(x, y, 'k')
            hold on
        end
    end
    
    hold on 
    set(gca, 'ydir', 'normal');
    hold on
    plot(xVortexCore(vtx), yVortexCore(vtx), 'wx', 'MarkerSize', 10, 'LineWidth', 1)
    axis equal
    axis([x_min x_max y_min y_max])
    cb = colorbar;
    set(get(cb, 'ylabel'), 'String', '\it \omega \rm (s^{-1})');
    
    % Set font and figure position
    set(gca, 'FontName', 'Times New Roman')
    set(gca, 'FontSize', 12)
    figs1 = figure(10 + vtx);
    style = hgexport('factorystyle');
    style.Bounds = 'tight';
    hgexport(figs1, '-clipboard', style, 'applystyle', true);
    drawnow;
    set(figure(10 + vtx), 'Position', [400 200 500 400])    
end
