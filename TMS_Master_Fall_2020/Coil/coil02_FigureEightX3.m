%   This script creates the mesh (both CAD surface mesh and a computational
%   wire grid) for a figure eight double planar coil with 1 A of total current
%   current. The output is saved in the binary file coil.mat and includes:
%   strcoil.Pwire(:, 3) - set of nodes for all wires 
%   strcoil.Ewire(:, 2) - set of edges or current dipoles for all wires
%   (current flows from the first edge node to the second edge node)
%   strcoil.Swire{:, 1} - current strength weight for every elementary
%   dipole asssuring that the total conductor current through any
%   cross-section is 1 A.

%   Copyright SNM 2018-2020

clear all; %#ok<CLALL>
if ~isunix
    s = pwd; addpath(strcat(s(1:end-5), '\Engine'));
else
    s = pwd; addpath(strcat(s(1:end-5), '/Engine'));
end

%   The coil is in the form of four interconnected
%   spiral arms. One vertical interconnection is ignored. The conductor
%   centerline model is given first
%   Two interconnected spiral arms
theta = [0*pi:pi/50:10*pi];
a0 = 0.01; b0 = 0.0010;
r = a0 + b0*theta;                          %   Archimedean spiral
x1 = r.*cos(theta);                         %   first half
y1 = r.*sin(theta);                         %   first half
x2 = 2*x1(end) - x1(end-1:-1:1);            %   second half
y2 = 2*y1(end) - y1(end-1:-1:1);            %   second half
X1 = [x1 x2]; Y1 = [y1 y2];                 %   join both halves
X1 = X1 - mean(X1);                         %   center the curve
L = length(X1)-1; arg = [1:L+1] -L/2 - 6;   
Z1 = 0.004*atan(arg);                       %   third argument z   
%   Two other interconnected spiral arms
X2 = - X1; Y2 = +Y1;                        %   use mirror refection here 
L = length(X1)-1; arg = [1:L+1] -L/2 + 6;   %   third argument z
Z2 = 0.004*atan(arg);
                
% plot3(X1, Y1, Z1, '*-r', X2, Y2, Z2, '*-b'); 
% axis equal; grid on; title('Conductor centerline')

%   Other parameters
a    = 4.0e-3;     %   z-side, m  (for a rectangle cross-section)
b    = 2.5e-3;     %   x-side, m  (for a rectangle cross-section)
M    = 16;         %   number of cross-section subdivisions 
flag = 2;          %   rectangular cross-section    
sk   = 1;          %   surface current distribution (skin layer)

%   Construct the mesh for first pair of interconnected spiral arms                      
Pcenter(:, 1) = X1';
Pcenter(:, 2) = Y1';
Pcenter(:, 3) = Z1;
%   Create CAD and wire models for the single conductor
strcoil1       = meshwire(Pcenter, a, b, M, flag, sk);
[P1, t1]       = meshsurface(Pcenter, a, b, M, flag);  %   CAD mesh (optional, slow)        

%   Construct the mesh for second pair of interconnected spiral arms
%   Construct the coil                       
Pcenter(:, 1) = X2';
Pcenter(:, 2) = Y2';
Pcenter(:, 3) = Z2;
%   Create CAD and wire models for the single conductor
strcoil2       = meshwire(Pcenter, a, b, M, flag, sk);
[P2, t2]       = meshsurface(Pcenter, a, b, M, flag);  %   CAD mesh (optional, slow)        

%   Combine both meshes together
strcoil.Ewire       = [strcoil1.Ewire; strcoil2.Ewire+size(strcoil1.Pwire, 1)];
strcoil.Pwire       = [strcoil1.Pwire; strcoil2.Pwire]; 
strcoil.Pwire(:, 3) = strcoil.Pwire(:, 3) - min(strcoil.Pwire(:, 3));
strcoil.Swire       = [strcoil1.Swire; strcoil2.Swire];
t   = [t1; t2+size(P1, 1)];
tind= [ones(size(t1, 1), 1); 2*ones(size(t2, 1), 1)]; 
P   = [P1; P2];
P(:, 3) = P(:, 3) - min(P(:, 3)); 
        
%   Display CAD and wire models
bemf1_graphics_coil_CAD(P, t, 0);
%bemf1_graphics_coil_wire(strcoil, [0 1 0]); 
view(-4, 24);
 
save('coil', 'strcoil');
save('coilCAD', 'P', 't', 'tind');  %   optional, slow
