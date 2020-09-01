function [ ] = bemf2_graphics_surf_field(P, t, FQ, Indicator, tissuenumber) 
%   Surface field graphics:  plot a field quantity FQ at the surface of a
%   brain compartment with the number "tissuenumber"
%
%   Copyright SNM 2017-2020

    ind = tissuenumber;                   
    t0  = t(Indicator==ind, :);   
    NumberOfTrianglesInShell = size(t0, 1);    
    patch('faces', t0, 'vertices', P, 'FaceVertexCData', FQ, 'FaceColor', 'flat', 'EdgeColor', 'none', 'FaceAlpha', 1.0);                   
    colormap jet; brighten(0.33); colorbar; 
    camlight; lighting phong;
    axis 'equal';  axis 'tight';      
    xlabel('x, m'); ylabel('y, m'); zlabel('z, m');
    set(gcf,'Color','White');    
end