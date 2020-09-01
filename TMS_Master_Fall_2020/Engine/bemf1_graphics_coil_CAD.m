function [ ] = bemf1_graphics_coil_CAD(P, t, flag) 
%   Coil 2D/3D plot with several options
%
%   Copyright SNM 2017-2020

    p = patch('vertices', P, 'faces', t);
    if flag == 0    %   non-transparent coil
        p.FaceColor = [1 0.75 0.65]; % [0.72 0.45 0.2];  
        p.EdgeColor = 'none';
        p.FaceAlpha = 1.0;
        daspect([1 1 1]);
        camlight('headlight'); lighting flat;   
    else
        p.FaceColor = [0.72 0.45 0.2];  
        p.EdgeColor = 'none';
        p.FaceAlpha = 0.20;      
    end
    xlabel('x, m'); ylabel('y, m'); zlabel('z, m');
    set(gcf,'Color','White');    
end