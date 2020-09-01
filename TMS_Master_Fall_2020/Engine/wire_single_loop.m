function [Pwire, Ewire] = wire_single_loop(Rx, Ry, a, b, N, M, type)
%   This function creates the structured wire grid for a single
%   ellipsoidal or circular loop with semi-axes Rx, Ry
%   a - cross-section conductor radius in m when option type='circle' is chosen for
%   the wire cross-section
%   a, b  - conductor rectangle sides in m when option type='rect' is chosen for
%   the wire cross-section
%   N - number of subdivisions along the loop circumference, 
%   M - number of cross-section subdivisions (approximate number for 
%   a rectangular cross-section) or the number of individual wires
%   The loop has the normal vector in the z-direction
%   The loop has M wires with N segments each
%   Pwire - array of output nodes
%   Ewire - array of output edges

%   Copyright SNM 2017-2020
%   The Athinoula A. Martinos Center for Biomedical Imaging, Massachusetts General
%   Hospital & ECE Dept., Worcester Polytechnic Inst.
    
    %   Step I
    if strcmp(type, 'circle')
        %   Create the structured grid for a cross-section in the form of a
        %   circle of radius a in the xz-plane with M subdivisions
        theta = [0:M-1]/M*(2*pi);
        xc = a*cos(theta);
        zc = a*sin(theta);    
    end
    if strcmp(type, 'rect')
        %   Create the structured grid for a cross-section in the form of a
        %   rectangle axb in the xz-plane with M subdivisions
        M4 = round(M/4);
        x = linspace(-a/2, +a/2, M4+1);    %   uniform grid
        z = linspace(-b/2, +b/2, M4+1);    %   uniform grid
        xc = [x, x(end)*ones(1, length(x)-2) x, x(1)*ones(1, length(x)-2)];
        zc = [z(1)*ones(1, length(x)), z(2:end-1), z(end)*ones(1, length(x)), z(2:end-1)];  
    end
    
    %   STEP II
    %   Create the complete array of M wires and associated edges
    S = N-1;
    Pwire = zeros(M*N, 3);
    Ewire = zeros(M*S, 2);
    phi   = linspace(0, 2*pi, N)';
    for m = 1:M     %   loop over wires
        Ewire(S*(m-1)+1:S*m, :) = [N*(m-1)+1:N*m-1; N*(m-1)+2:N*m]'; 
        RadX    = Rx + xc(m);
        RadY    = Ry + xc(m);
        for n = 1:N %   loop over points in a single wire
            Pwire(n+N*(m-1), 1) = RadX*cos(phi(n));
            Pwire(n+N*(m-1), 2) = RadY*sin(phi(n));
            Pwire(n+N*(m-1), 3) = zc(m);
        end
    end     
end
    