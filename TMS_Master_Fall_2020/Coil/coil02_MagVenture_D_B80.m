%   This script creates the mesh (both CAD surface mesh and a computational
%   wire grid) for a MagVenture figure eight bent coil with 1 A of total
%   current
%   The output is saved in the binary file coil.mat and includes:
%   strcoil.Pwire(:, 3) - set of nodes for all wires 
%   strcoil.Ewire(:, 2) - set of edges or current dipoles for all wires
%   (current flows from the first edge node to the second edge node)
%   strcoil.Swire{:, 1} - current strength weight for every elementary
%   dipole asssuring that the total conductor current through any
%   cross-section is 1 A.
%
%   Copyright SNM 2020

clear all; %#ok<CLALL>
if ~isunix
    s = pwd; addpath(strcat(s(1:end-5), '\Engine'));
else
    s = pwd; addpath(strcat(s(1:end-5), '/Engine'));
end

%   The coil is in the form of two interconnected spiral arms. The
%   conductor centerline model is given first

%%   First arm
theta = [8*pi/2:pi/50:12*pi];
a0 = 0.024; b0 = 0.00061;
r = a0 + b0*theta;                  %   Archimedean spiral
x1 = r.*cos(theta);                 %   first half
y1 = r.*sin(theta);                 %   first half
x2 = 2*x1(end) - x1(end-1:-1:1);    %   second half
y2 = 2*y1(end) - y1(end-1:-1:1);    %   second half
x = [x1 x2];  y = [y1 y2];          %   join both halves
x = x - mean(x);                    %   center the curve

% plot(x, y, '*-'); axis equal; grid on; title('Conductor centerline')

%   Other parameters
a    = 6.0e-3;     %   z-side, m  (for a rectangle cross-section)
b    = 2.0e-3;     %   x-side, m  (for a rectangle cross-section)
M    = 20;          %   number of cross-section subdivisions 
flag = 2;           %   rect. cross-section    
sk   = 1;           %   surface current distribution (skin layer)

%   Create CAD and wire models for the single conductor
Pcenter(:, 1) = x';
Pcenter(:, 2) = y';
Pcenter(:, 3) = -a/2;
strcoil1       = meshwire(Pcenter, a, b, M, flag, sk);
[P1, t1]        = meshsurface(Pcenter, a, b, M, flag);  %   CAD mesh (optional, slow)     
tind1          = 1*ones(size(t1, 1), 1);

%%   Second arm
theta = [8*pi/2:pi/50:10*pi];
a0 = 0.027; b0 = 0.00061;
r = a0 + b0*theta;                  %   Archimedean spiral
x1 = r.*cos(theta);                 %   first half
y1 = r.*sin(theta);                 %   first half
x2 = 2*x1(end) - x1(end-1:-1:1);    %   second half
y2 = 2*y1(end) - y1(end-1:-1:1);    %   second half
x = [x1 x2];  y = [y1 y2];          %   join both halves
x = x - mean(x);                    %   center the curve

% plot(x, y, '*-'); axis equal; grid on; title('Conductor centerline')

%   Other parameters
a    = 6.0e-3;     %   z-side, m  (for a rectangle cross-section)
b    = 2.0e-3;     %   x-side, m  (for a rectangle cross-section)
M    = 20;          %   number of cross-section subdivisions 
flag = 2;           %   rect. cross-section    
sk   = 1;           %   surface current distribution (skin layer)

%   Create CAD and wire models for the single conductor
clear Pcenter;
Pcenter(:, 1) = x';
Pcenter(:, 2) = -y';
Pcenter(:, 3) = +a/2;
strcoil2       = meshwire(Pcenter, a, b, M, flag, sk);
[P2, t2]      = meshsurface(Pcenter, a, b, M, flag);  %   CAD mesh (optional, slow)     
tind2         = 2*ones(size(t2, 1), 1);

%%  Create the complete model
Ewire       = [];
Pwire       = []; 
Swire       = [];
Pa          = [];
ta          = []; 

%   Construct two CAD and wire models 
strcoil.Swire       = [strcoil1.Swire; -strcoil2.Swire];
strcoil.Ewire       = [strcoil1.Ewire; strcoil2.Ewire+size(strcoil1.Pwire, 1)];
strcoil.Pwire       = [strcoil1.Pwire; strcoil2.Pwire]; 
strcoil.Pwire(:, 3) = strcoil.Pwire(:, 3) - min(strcoil.Pwire(:, 3));

t          = [t1; t2+size(P1, 1)];
tind       = [tind1; tind2];
P          = [P1; P2];
P(:, 3)    = P(:, 3) - min(P(:, 3)); 

%   Deform the entire structure
alpha = pi/6;
strcoil.Pwire(:, 3) = strcoil.Pwire(:, 3) - sin(alpha)*abs(strcoil.Pwire(:, 1));
strcoil.Pwire(:, 1) =                       cos(alpha)*strcoil.Pwire(:, 1);
P(:, 3)             = P(:, 3) - sin(alpha)*abs(P(:, 1));
P(:, 1)             =           cos(alpha)*P(:, 1);

%   Display CAD and wire models
bemf1_graphics_coil_CAD(P, t, 0);
%bemf1_graphics_coil_wire(strcoil, [0 1 0]); 
view(20, 20);
 
save('coil', 'strcoil');
save('coilCAD', 'P', 't', 'tind');  %   optional, slow