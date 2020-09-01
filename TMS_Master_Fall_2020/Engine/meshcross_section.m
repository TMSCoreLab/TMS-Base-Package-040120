function [P, e] = meshcross_section(a, b, normal, M, flag)
    %   Creates a structured edge grid P, e for a perimeter 
    %   of a base ellipse (flag = 1) or rectangle (flag = 2) with
    %   - major axis/side a (say long one; always in the z direction);
    %   - minor axis/side b (say short one; originally in the x direction); 
    %   - normal vector of the ellipse given by unit normal = [nx, ny, 0];
    %   and (approximately for rectangle) M edges. 
    %   The grid is centered at the origin.    
    %   Copyright SNM 2018-2020
    
    if flag == 1
        x = b/2*cos(2*pi*[0:M-1]/M);
        z = a/2*sin(2*pi*[0:M-1]/M);
    else
        M4 = round(M/4);
        x = linspace(-b/2, +b/2, M4);    %   uniform grid
        z = linspace(-a/2, +a/2, M4);    %   uniform grid
        xc = [x, ...
              x(end)*ones(1, length(x)-2),...
              x(end:-1:1),...
              x(1)*ones(1, length(x)-2)];
        zc = [z(1)*ones(1, length(x)),...
              z(2:end-1),...
              z(end)*ones(1, length(x)),...
              z(end-1:-1:2)];  
        x = xc; z = zc;  
    end
        
    P(:, 1) = x';
    P(:, 2) = 0;
    P(:, 3) = z';
    e(:, 1) = [1:size(P, 1)  ]';
    e(:, 2) = [2:size(P, 1) 1]';
    
    %   Rotate mesh
    angle    =  acos(normal(2)/norm(normal));
    if normal(1)>0; angle = 2*pi - angle; end;
    P        =  meshrotate2(P, [0 0 1], angle);
end
    