function [P, t, nodesadded] = meshfill(p, normal)
    %   Creates an inner mesh dt for an arbitrarily oriented planar polygon
    %   Inputs:
    %   p       - polygon in 3D, [:,3]
    %   normal  - normal vector to the polygon surface
    %   Outputs:
    %   dt          - output mesh
    %   nodesadded  - number of added nodes      
    %   Copyright SNM 2018-2020
    
    %  Internal parameter - number of mesh refinement steps
    N = 2;      
    p                   = [p; mean(p, 1)];      %   add center node
    nodesadded          = 1;    
    [dummy, index]      = max(abs(normal));             
    for m = 1:N
        if index ==1  P = p(:, [2 3]); end;    %   to do 2D Delaunay
        if index ==2  P = p(:, [1 3]); end;    %   to do 2D Delaunay
        if index ==3  P = p(:, [1 2]); end;    %   to do 2D Delaunay
        dt               = delaunayTriangulation(P);  %   2D Delaunay
        if m<N 
            C = meshtricenter(p, dt.ConnectivityList);
            p = [p; C];
            nodesadded = nodesadded + size(C, 1);
        end
    end
    P = p;
    t = dt.ConnectivityList;
end
    