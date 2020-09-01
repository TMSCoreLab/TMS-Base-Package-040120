%   This script creates the mesh (both CAD surface mesh and a computational
%   wire grid) for a MagVenture figure eight planar coil with 1 A of total
%   current
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

%   The coil is in the form of two non-interconnected spiral arms. The
%   conductor centerline model is given first
turns = 32;
theta = [0:pi/20:turns*2*pi-pi/2];
a0 = 0.0115; b0 = (23e-3 - a0)/(turns*2*pi-pi/2);
r = a0 + b0*theta;                 %   Archimedean spiral
x = r.*cos(theta);                 %   first half
y = r.*sin(theta);                 %   first half

% plot(x, y, '*-'); axis equal; grid on; title('Conductor centerline');
% return

%   Other parameters
a    = 15e-3;       %   z-side, m  (for a rectangle cross-section)
b    = 0.2e-3;      %   x-side, m  (for a rectangle cross-section)
M    = 20;          %   number of cross-section subdivisions 
flag = 2;           %   rect. cross-section    
sk   = 0;           %   Litz wire

%   Create CAD and wire models for the single conductor
Pcenter(:, 1) = x';
Pcenter(:, 2) = y';
Pcenter(:, 3) = a/2;
strcoil       = meshwire(Pcenter, a, b, M, flag, sk);
[P, t]        = meshsurface(Pcenter, a, b, M, flag);  %   CAD mesh (optional, slow)     
tind          = 1*ones(size(t, 1), 1);

Ewire       = [];
Pwire       = []; 
Swire       = [];
Pa          = [];
ta          = []; 

%   Construct two CAD and wire models 
strcoil.Swire       = [strcoil.Swire; strcoil.Swire];
strcoil.Ewire       = [strcoil.Ewire; strcoil.Ewire+size(strcoil.Pwire, 1)];
Pwire1              = strcoil.Pwire;
Pwire2              = strcoil.Pwire;
Pwire1(:, 1)        = -Pwire1(:, 1);
Pwire1(:, 1)        = Pwire1(:, 1) - 23e-3;
Pwire2(:, 1)        = Pwire2(:, 1) + 23e-3;
strcoil.Pwire       = [Pwire1; Pwire2]; 
strcoil.Pwire(:, 3) = strcoil.Pwire(:, 3) - min(strcoil.Pwire(:, 3));

t          = [t; t+size(P, 1)];
tind       = [tind; 2*tind];
P1         = P;
P2         = P;
P1(:, 1)   = -P1(:, 1);
P1(:, 1)   = P1(:, 1) - 23e-3;
P2(:, 1)   = P2(:, 1) + 23e-3;
P          = [P1; P2]; 
P(:, 3)    = P(:, 3) - min(P(:, 3)); 

%   Display CAD and wire models
bemf1_graphics_coil_CAD(P, t, 0);
%bemf1_graphics_coil_wire(strcoil, [0 1 0]); 
view(5, 65);
 
save('coil', 'strcoil');
save('coilCAD', 'P', 't', 'tind');  %   optional, slow