function [P, t] = meshsurface(Pcenter, a, b, M, flag)
%   Outputs a 2-manifold P, t mesh for a single arbitrarily 
%   bent conductor. The conductor could be either open or closed. In the
%   last case, the start point and the end point must coincide.
%   Inputs:
%   Pcenter - centerline of the conductor in 3D [:, 3];
%   a - major axis/side (always in the z direction) of conductor cross-section;
%   b - minor axis/side (always in the z direction) of conductor cross-section;
%   M - number of cross-section subdivisions (approximate for rectangular
%   cross-section);
%   flag is equal to one for the elliptical cross-section and equals two for
%   the rectangular cross-section
%   Outputs:
%   P  -  P-aray of surface mesh vertices
%   t  -  t-array of surface triangular facets
%   SNM 2018-2020
    
    %%  Create 2-manifold surface mesh
    Nodes           = size(Pcenter, 1); 
    Closed          = norm(Pcenter(1, :)-Pcenter(end, :))<1024*eps;
    PathVector      = Pcenter(2:end, :) - Pcenter(1:end-1, :);       
    %   Add nodes/triangles for conductor side surface   
    t           = [];    
    for m = 1:Nodes - 1        
        if m == 1 % bottom only
            UnitPathVector  = PathVector(1, :) + Closed*PathVector(end, :);
            UnitPathVector  = UnitPathVector/norm(UnitPathVector);                  
            [pbottom, e]    = meshcross_section(a, b, UnitPathVector, M, flag);
            NE              = size(e, 1);   % number of nodes/edges in the cross-section, global for the entire code
            pbottom         = pbottom + repmat(Pcenter(m, :), NE, 1);                 
            P               = pbottom;     
        end    
        if m >= 1 & m < Nodes - 1
            UnitPathVector  = PathVector(m, :) + PathVector(m+1, :);
            UnitPathVector  = UnitPathVector/norm(UnitPathVector);  
            [ptop, e]       = meshcross_section(a, b, UnitPathVector, M, flag);  
            ptop            = ptop + repmat(Pcenter(m+1, :), NE, 1); 
            P               = [P; ptop];                                
        end  
        if m == Nodes - 1 % top only
            UnitPathVector  = PathVector(end, :) + Closed*PathVector(1, :);
            UnitPathVector  = UnitPathVector/norm(UnitPathVector);   
            [ptop, e]   = meshcross_section(a, b, UnitPathVector, M, flag);
            ptop        = ptop + repmat( Pcenter(m+1, :), NE, 1);  
            P           = [P; ptop];      
        end        
        %   Local connectivity: bottom to top
        t1(:, 1:2)  = e;                        %   Lower nodes        
        t1(:, 3)    = e(:, 1) + NE;             %   Upper node
        t2(:, 2:-1:1)  = e       + NE;             %   Upper nodes
        t2(:, 3)    = e(:, 2);                  %   Lower node
        ttemp       = [t1; t2];
        t           = [t; ttemp+NE*(m-1)];
    end

    %%   Add caps for a non-closed conductor
    if ~Closed  
        %   Add nodes/triangles for the start cap        
        PathVector1              = Pcenter(2, :) - Pcenter(1, :);
        [p1, tcapstart, added]  = meshfill(pbottom, PathVector1);         
        t1                      = tcapstart;
        NodesSides              = size(P, 1)-NE;
        for m = 1:size(tcapstart, 1)
            index               = find(tcapstart(m, :)>NE);
            for n = 1:size(index, 2)
                tcapstart(m, index(n)) = tcapstart(m, index(n)) + NodesSides;    
            end
        end
        P           = [P; p1(NE+1:end, :)];
        t           = [tcapstart; t];  
        %   Add nodes/triangles for the end cap
        [p2, tcapend, ~]    = meshfill(ptop, UnitPathVector);
        tcapend0            = tcapend;
        NodesSides1         = size(P, 1)-NE-added;
        NodesSides2         = size(P, 1)-NE;
        for m = 1:size(tcapend, 1)
            index               = find(tcapend0(m, :)<=NE);
            for n = 1:size(index, 2)
                tcapend(m, index(n)) = tcapend(m, index(n)) + NodesSides1;    
            end
            index               = find(tcapend0(m, :)>NE);
            for n = 1:size(index, 2)
                tcapend(m, index(n)) = tcapend(m, index(n)) + NodesSides2;    
            end    
        end
        P        = [P; p2(NE+1:end, :)];
        t        = [t; tcapend];    
    end
    %   Condition the final surface mesh
    [P, t] = fixmesh(P, t);
end

