function [Pwire, Ewire, Swire, P, t, tind] = meshcoil(x0, y0, z0, M, N, a, b, flag, sk)
%   Create the mesh (both CAD surface mesh and a computational wire grid)
%   for a simple z-oriented coil with elliptical/spherical turns
%   characterized by centerline intersections x0, y0, z0 with xz and yz
%   planes and 1 A of total current per conductor.
%   Inputs:
%   x0, y0, z0 - row vectors which are centerline intersections with xz and yz
%   planes. Number of turns is equal to the length of x0, y0, z0
%   M - number of cross-section perimeter subdivisions (approximate for
%   rectangular cross-section)
%   N - number of loop perimeter subdivisions
%   a - major axis/side (always in the z direction) of conductor cross-section
%   b - minor axis/side (always in the z direction) of conductor cross-section
%   flag is equal to one for the elliptical cross-section and equals two for
%   the rectangular cross-section
%   parameter sk equals to zero for uniform current distribution (Litz wire)
%   or to 1 for the skin effect (bulk of the current flows close to the
%   surface)
%   Outputs:
%   W.Pwire - nodes of elementary wires inside the conductor
%   W.Ewire - edges of elementary wires inside the conductor
%   W.Swire - weights of elementary wire segments given total current of 1A
%   P  -  P-aray of surface mesh vertices
%   t  -  t-array of surface triangular facets

%   Copyright SNM 2018-2020

    %   Construct the coil (centered)
    Pwire = []; Ewire = []; Swire = [];     %   computational wire grid 
    P = []; t = [];                         %   CAD mesh (optional, slow)
    tind = [];
    for m = 1:length(x0)
        %   Define the loop centerline        
        x = x0(m)*cos(2*pi*[0:N]/N);
        y = y0(m)*sin(2*pi*[0:N]/N);
        Pcenter(:, 1) = x';
        Pcenter(:, 2) = y';
        Pcenter(:, 3) = z0(m);
        %   Define the loop cross-section  
        %   Create CAD and wire models
        W               = meshwire(Pcenter, a, b, M, flag, sk);
        [P1, t1]        = meshsurface(Pcenter, a, b, M, flag);            
        Ewire = [Ewire; W.Ewire+size(Pwire, 1)];
        Pwire = [Pwire; W.Pwire];    
        Swire = [Swire; W.Swire];
        t     = [t; t1+size(P, 1)];
        P     = [P; P1];
        tind  = [tind; m*ones(size(t1, 1), 1)]; 
    end
    [P, t] = fixmesh(P, t);   
end

