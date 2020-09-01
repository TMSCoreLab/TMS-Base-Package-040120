function W = meshwire(Pcenter, a, b, M, flag, sk)
%   Outputs structure W with the equivalent FMM computational wire grid
%   Inputs:
%   Pcenter - centerline of the conductor in 3D [:, 3];
%   a - major axis/side (always in the z direction) of conductor cross-section;
%   b - minor axis/side (always in the z direction) of conductor cross-section;
%   M - number of cross-section subdivisions (approximate for rectangular
%   cross-section);
%   flag is equal to one for the elliptical cross-section and equals two for
%   the rectangular cross-section
%   parameter sk equals to zero for uniform current distribution (Litz wire)
%   or to 1 for the skin effect (bulk of the current flows close to the
%   surface)
%   Outputs:
%   W.Pwire - nodes of elementary wires inside the conductor
%   W.Ewire - edges of elementary wires inside the conductor
%   W.Swire - weights of elementary wire segments given total current of 1A
%
%   SNM 2018-2020
    
    %%  Create computational wire grid Pwire, Ewire, Swire based on the mesh for the starting cap 
    Nodes           = size(Pcenter, 1);
    Closed          = norm(Pcenter(1, :)-Pcenter(end, :))<1024*eps;
    PathVector      = Pcenter(2:end, :) - Pcenter(1:end-1, :);       

    %   Create triangular mesh for the start cap 
    UnitPathVector  = PathVector(1, :) + Closed*PathVector(end, :);
    UnitPathVector  = UnitPathVector/norm(UnitPathVector);                  
    [pbottom, e]    = meshcross_section(a, b, UnitPathVector, M, flag); 
    NE              = size(e, 1);   % number of nodes/edges in the cross-section, global for the entire code
    pbottom         = pbottom + repmat(Pcenter(1, :), NE, 1);    
    [p, t, ~]       = meshfill(pbottom, PathVector(1, :));             

    %   Process triangular mesh for the start cap
    edges   = meshconnee(t); 
    areas   = meshareas(p, t);    
    points  = meshtricenter(p, t);      
    %   Determine weights including border triangles (if necessary)
    if sk == 0
       weights     = areas/sum(areas);  % uniform curent distribution/weights
    else
        at    = meshconnet(t, edges, "nonmanifold"); % triangles attached to every edge
        tborder = [];
        for m = 1:length(at)
            if at{m}(2)==0
                tborder = [tborder at{m}(1)];
            end
        end
        weights  = areas(tborder)/sum(areas(tborder));
        points   = points(tborder, :); 
    end
    
    %   Allocate and fill out wire arays
    wires       = size(points, 1);
    Pwire       = zeros((Nodes-0)*wires, 3);
    Ewire       = zeros((Nodes-1)*wires, 2);
    Swire       = repmat(weights, Nodes-1, 1);    
    arg         = zeros(Nodes, wires);
    for m = 1:Nodes  
        arg(m, :)  = [1:wires]+(m-1)*wires;
        if m == 1 %   define the initial set of nodes
            UnitPathVector  = PathVector(1, :) + Closed*PathVector(end, :);
            UnitPathVector  = UnitPathVector/norm(UnitPathVector);
            angle0          = acos(UnitPathVector(2));           
            if UnitPathVector(1)>0; angle0 = 2*pi - angle0; end;              
            base            = points - repmat(Pcenter(m, :), wires, 1);
        end    
        if m > 1   %   define the following sets of nodes
            if m < Nodes
                UnitPathVector  = PathVector(m, :) + PathVector(m-1, :);
            else 
                UnitPathVector  = PathVector(end, :) + Closed*PathVector(1, :);
            end
            UnitPathVector  = UnitPathVector/norm(UnitPathVector);  
            angle1          = acos(UnitPathVector(2));          %   restore the rotation 
            if UnitPathVector(1)>0; angle1 = 2*pi - angle1; end; 
            points          =  meshrotate2(base, [0 0 1], angle1-angle0); 
            points          =  points + repmat(Pcenter(m, :), wires, 1);
            Ewire(arg(m-1, :), :)   = [arg(m-1, :); arg(m, :)]';       
        end 
        Pwire(arg(m, :), :)     = points;
    end
    W.Pwire = Pwire; W.Ewire = Ewire; W.Swire = Swire;
end

