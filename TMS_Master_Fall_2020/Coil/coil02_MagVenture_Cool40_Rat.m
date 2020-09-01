%   This script creates the mesh (both CAD surface mesh and a computational
%   wire grid) for a MagVenture Cool40_Rat coil with 1 A of total current
%   The output is saved in the binary file coil.mat and includes:
%   strcoil.Pwire(:, 3) - set of nodes for all wires strcoil.Ewire(:, 2) -
%   set of edges or current dipoles for all wires (current flows from the
%   first edge node to the second edge node) strcoil.Swire{:, 1} - current
%   strength weight for every elementary dipole asssuring that the total
%   conductor current through any cross-section is 1 A.
%
%   Copyright SNM 2018-2020

clear all; %#ok<CLALL>
if ~isunix
    s = pwd; addpath(strcat(s(1:end-5), '\Engine'));
else
    s = pwd; addpath(strcat(s(1:end-5), '/Engine'));
end

strcoil.Swire       = [];
strcoil.Ewire       = [];
strcoil.Pwire       = [];
P                   = [];
t                   = [];
tind                = [];
%   Creating the coil turns in a loop
for m  = 1:36
    offset = 0;
    if m>12
        offset = -4.3e-3;
        m = m -12;
    end
    if m>24
        offset = -2*4.3e-3;
        m = m -24;
    end
    %   The coil is in the form of many deformed ellipsoidal loops. The
    %   conductor centerline model is given first
    delta0  = 0.065;
    delta1  = 0.065;
    a0      = 19.0e-3*(1-delta0*(m-1))-0.5e-3;     %   minor radius xy plane
    a1      = 22.5e-3*(1-delta1*(m-1))-0.5e-3;     %   major radius xy plane
    R       = 34e-3;                        %   bending radius yz plane
    theta   = linspace(0, 2*pi, 100);
    y = a1.*cos(theta);                 %   planar ellipse
    x = a0.*sin(theta);                 %   planar ellipse
    z = sqrt(R^2 - x.^2) - R + offset;  %   deformed ellipse xz plane
    %   Other parameters
    a    = 3e-3;     %   z-side, m  (for a rectangle cross-section)
    b    = 0.5e-3;     %   x-side, m  (for a rectangle cross-section)
    M    = 20;          %   number of cross-section subdivisions 
    flag = 2;           %   rect. cross-section    
    sk   = 1;           %   surface current distribution (skin layer)

    %   Create CAD and wire models for the single conductor
    Pcenter(:, 1) = x';
    Pcenter(:, 2) = y';
    Pcenter(:, 3) = z';
    strcoil_temp    = meshwire(Pcenter, a, b, M, flag, sk);
    [Ptemp, ttemp]  = meshsurface(Pcenter, a, b, M, flag);  %   CAD mesh (optional, slow)     
    tindtemp       = 1*ones(size(t, 1), 1);

    %   Accumulate loops 
    strcoil.Swire       = [strcoil.Swire; strcoil_temp.Swire];
    strcoil.Ewire       = [strcoil.Ewire; strcoil_temp.Ewire+size(strcoil.Pwire, 1)];
    strcoil.Pwire       = [strcoil.Pwire; strcoil_temp.Pwire];     
    
    t          = [t; ttemp+size(P, 1)];    
    tind       = [tind; tindtemp];    
    P          = [P; Ptemp];  
end

%   Display CAD and wire models
bemf1_graphics_coil_CAD(P, t, 0);
%bemf1_graphics_coil_wire(strcoil, [0 1 0]); 
view(-4, 24);
 
save('coil', 'strcoil');
save('coilCAD', 'P', 't', 'tind');  %   optional, slow