function [figs3, figs4] = plotProfiles(numROIs, velocityData, xVortexCore, yVortexCore, radiiVortexCore, drift_vy, rotationSign)
% PLOTPROFILES - Visualize velocity and vorticity profiles for a line crossing the core of each vortex.
%
% Description:
%   This function generates figures to visualize velocity and vorticity profiles along a line
%   that crosses the core of each vortex. It also compares the sample vortex with ideal Rankine
%   and Burgers vortices.
%
% Inputs:
%   - numROIs: Number of vortex cores (integer).
%   - velocityData: Structure containing velocity field data.
%   - xVortexCore: x-coordinate of vortex cores.
%   - yVortexCore: y-coordinate of vortex cores.
%   - radiiVortexCore: Radius of vortex cores.
%   - drift_vy: Drift velocity in the y-direction for each vortex core.
%   - rotationSign: Array indicating the rotation sign for each vortex core.
%
% Outputs:
%   - figs3: Handle to the figure displaying vorticity profiles.
%   - figs4: Handle to the figure displaying velocity profiles.
%
%--------------------------------------------------------------------------

for vtx = 1:numROIs
    % Initialize variables
    clear Box_L
    clear theta_wkind
    clear vb
    clear vr
    clear wb
    clear wr 
    
    % Define constants
    box_fact = 1.65;
    gridSpacing = velocityData.x(end) - velocityData.x(end - 1);
    
    % Compute boundaries of the vortex along CLSc
    Box_L = floor(box_fact * radiiVortexCore(vtx) / gridSpacing);
    Nn = 2 * Box_L;
    
    [~, y_indexc] = min(abs(velocityData.y - yVortexCore(vtx))); 
    [~, x_indexc] = min(abs(velocityData.x - xVortexCore(vtx))); 

    Lbi = x_indexc - Box_L;
    Rbi = x_indexc + Box_L;

    if Rbi > size(velocityData.x, 2)
        Rbi = size(velocityData.x, 2);
        Nn = (floor((Rbi - Lbi) / 2)) * 2;
    end
    if Lbi < 1
        Lbi = 1;
        Nn = Rbi - Lbi + 1;
    end

    % Find v.vy and w values along CLSc
    v_CSLc{vtx} = velocityData.vy(y_indexc, Lbi:Rbi) - drift_vy(vtx); 
    w_CSLc{vtx} = velocityData.rot(y_indexc, Lbi:Rbi);

    box_fact2 = 0.75;
    
    % Build ideal Rankine and Burgers vortices
    vb = vortex(Nn + 1, radiiVortexCore(vtx) * 0.77, max(abs(w_CSLc{vtx})), 'burgers');
    vr = vortex(Nn + 1, radiiVortexCore(vtx) * 0.77, max(abs(w_CSLc{vtx})), 'rankine');
    wb = vec2scal(vb, 'curl');
    wr = vec2scal(vr, 'curl');
    
    wb.w = wb.w';
    vb.vy = vb.vy';
    vb.vx = vb.vx';
    wr.w = wr.w';
    vr.vy = vr.vy';
    vr.vx = vr.vx';

    vb_CSLc{vtx} = rotationSign{vtx} .* -vb.vy(round(Nn / 2), :);
    wb_CSLc{vtx} = rotationSign{vtx} .* wb.w(round(Nn / 2), :);
    
    % Dimensionless
    vb_CSLcd{vtx} = vb_CSLc{vtx} / max(abs(vb_CSLc{vtx}));
    wb_CSLcd{vtx} = wb_CSLc{vtx} / max(abs(wb_CSLc{vtx}));

    vr_CSLc{vtx} = rotationSign{vtx} .* -vr.vy(round(Nn / 2), :);
    wr_CSLc{vtx} = rotationSign{vtx} .* wr.w(round(Nn / 2), :);
    
    % Dimensionless
    vr_CSLcd{vtx} = vr_CSLc{vtx} / max(abs(vr_CSLc{vtx}));
    wr_CSLcd{vtx} = wr_CSLc{vtx} / max(abs(wr_CSLc{vtx}));
    
    % Dimensionless x
    dimlessx{vtx} = (velocityData.x(Lbi:Rbi) - xVortexCore(vtx)) / (radiiVortexCore(vtx));

    % Create figure for vorticity profiles
    figs3 = figure(30 + vtx);
    plot(dimlessx{vtx}, w_CSLc{vtx});
    hold on;
    plot(dimlessx{vtx}, wb_CSLc{vtx}, 'r--');
    plot(dimlessx{vtx}, wr_CSLc{vtx}, 'k-.');
    xlabel('\it (y-y_{c}) / r_{v}');
    ylabel('\it \omega \rm (s^{-1})');
    legend('Sample vortex', 'Burgers vortex', 'Rankine vortex');
    hold on;

    % Create figure for velocity profiles
    figs4 = figure(40 + vtx);
    plot(dimlessx{vtx}, v_CSLc{vtx});
    hold on;
    plot(dimlessx{vtx}, vb_CSLc{vtx}, 'r--');
    plot(dimlessx{vtx}, vr_CSLc{vtx}, 'k-.');
    xlabel('\it (y-y_{c}) / r_{v}');
    ylabel('\it v_{\it\theta} \rm (mm s^{-1})');
    legend('Sample vortex', 'Burgers vortex', 'Rankine vortex');
    hold on;

    % Set font and figure position
    set(gca, 'FontName', 'Times New Roman');
    style = hgexport('factorystyle');
    style.Bounds = 'tight';
    drawnow;
end

