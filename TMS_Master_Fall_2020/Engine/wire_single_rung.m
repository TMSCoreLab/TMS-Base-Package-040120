function [Pwire, Ewire] = wire_single_rung(L, a, b, N, M, type)
%   This function creates the structured centered wire grid for a single
%   straight rung with
%   a - cross-section radius in m when option type='circle' is chosen for
%   the rung cross-section
%   a, b  - rectangle sides in m when option type='rect' is chosen for
%   the rung cross-section
%   N - number of subdivisions along the rung length, 
%   M - number of cross-section subdivisions (approximate number for 
%   a rectangular cross-section)
%   The rung is directed along the z-axis
%   The rung has M wires with N segments each
%   Pwire - array of output nodes
%   Ewire - array of output edges

%   Copyright SNM 2017-2020
%   The Athinoula A. Martinos Center for Biomedical Imaging, Massachusetts General
%   Hospital & ECE Dept., Worcester Polytechnic Inst.
    
    %   Step I
    if strcmp(type, 'circle')
        %   Create the structured grid for a cross-section in the form of a
        %   circle of radius a in the xy-plane with M subdivisions
        theta = [0:M-1]/M*(2*pi);
        xc = a*cos(theta);
        yc = a*sin(theta);    
    end
    if strcmp(type, 'rect')
        %   Create the structured grid for a cross-section in the form of a
        %   rectangle axb in the xy-plane with M subdivisions
        M4 = round(M/4);
        x = linspace(-a/2, +a/2, M4+1);    %   uniform grid
        y = linspace(-b/2, +b/2, M4+1);    %   uniform grid
        xc = [x, x(end)*ones(1, length(x)-2) x, x(1)*ones(1, length(x)-2)];
        yc = [y(1)*ones(1, length(x)), y(2:end-1), y(end)*ones(1, length(x)), y(2:end-1)];  
    end
    
    %   STEP II 
    %   Create the complete array of M wires and associated edges
    S = N-1;
    Pwire = zeros(M*N, 3);
    Ewire = zeros(M*S, 2);
    z   = linspace(-L/2, +L/2, N)';
    for m = 1:M     %   loop over wires
        Ewire(S*(m-1)+1:S*m, :) = [N*(m-1)+1:N*m-1; N*(m-1)+2:N*m]'; 
        for n = 1:N %   loop over points in a single wire
            Pwire(n+N*(m-1), 1) = xc(m);
            Pwire(n+N*(m-1), 2) = yc(m);
            Pwire(n+N*(m-1), 3) = z(n);
        end
    end    
end
    