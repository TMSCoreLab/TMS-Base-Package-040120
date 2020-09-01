function [ ] = bemf1_graphics_coil_wire(strcoil, color) 
%   Coil 2D/3D plot with several options
%   type  = 0 for a 3D plot
%   type  = 1 for a planar plot in the coronal plane (XZ)
%   type  = 2 for a planar plot in the sagittal plane(YZ)
%   type  = 3 for a planar plot in the transverse plane(XY)
%   color = 'g':  color for wires
%
%   SNM 2017-2020

    lw = 2;
    patch('faces', strcoil.Ewire, 'vertices', strcoil.Pwire, 'EdgeColor', 'k', 'FaceColor', color, 'LineWidth', lw);
    axis 'equal';  axis 'tight';      
    xlabel('x, m'); ylabel('y, m'); zlabel('z, m');
    set(gcf,'Color','White');
end