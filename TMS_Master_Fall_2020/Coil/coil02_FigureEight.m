%   This script creates the mesh (both CAD surface mesh and a computational
%   wire grid) for a figure eight planar coil with 1 A of total current
%   The output is saved in the binary file coil.mat and includes:
%   strcoil.Pwire(:, 3) - set of nodes for all wires 
%   strcoil.Ewire(:, 2) - set of edges or current dipoles for all wires
%   (current flows from the first edge node to the second edge node)
%   strcoil.Swire{:, 1} - current strength weight for every elementary
%   dipole asssuring that the total conductor current through any
%   cross-section is 1 A.
%
%   Copyright SNM 2018-2020

clear all; %#ok<CLALL>
if ~isunix
    s = pwd; addpath(strcat(s(1:end-5), '\Engine'));
else
    s = pwd; addpath(strcat(s(1:end-5), '/Engine'));
end

%   The coil includes one single conductor in the form of two
%   interconnected spiral arms. The conductor centerline model is given
%   first
theta = [0:pi/25:10*pi];
a0 = 0.03; b0 = 0.0006;
r = a0 + b0*theta;                  %   Archimedean spiral
x1 = r.*cos(theta);                 %   first half
y1 = r.*sin(theta);                 %   first half
x2 = 2*x1(end) - x1(end-1:-1:1);    %   second half
y2 = 2*y1(end) - y1(end-1:-1:1);    %   second half
x = [x1 x2];  y = [y1 y2];          %   join both halves
x = x - mean(x);                    %   center the curve

%   Other parameters
a    = 0.002;       %   conductor diameter
M    = 16;          %   number of cross-section subdivisions 
flag = 1;           %   circular cross-section    
sk   = 1;           %   surface current distribution (skin layer)

%   Create CAD and wire models for the single conductor
Pcenter(:, 1) = x';
Pcenter(:, 2) = y';
Pcenter(:, 3) = a/2;
strcoil       = meshwire(Pcenter, a, a, M, flag, sk); % Wire model 
[P, t]        = meshsurface(Pcenter, a, a, M, flag);  % CAD mesh (optional, slow)     
tind          = ones(size(t, 1), 1); 

%   Display CAD and wire models
bemf1_graphics_coil_CAD(P, t, 0);
%bemf1_graphics_coil_wire(strcoil, [0 1 0]); 
view(-4, 24);
 
save('coil', 'strcoil');
save('coilCAD', 'P', 't', 'tind');  %   optional, slow
