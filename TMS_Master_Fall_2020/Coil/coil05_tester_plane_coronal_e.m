%   This script accurately computes and displays electric field sampled
%   over a surface (coronal plane) due to a coil via the plain FMM method
%
%   Copyright SNM 2018-2020

clear all; %#ok<CLALL>
if ~isunix
    s = pwd; addpath(strcat(s(1:end-5), '\Engine'));
else
    s = pwd; addpath(strcat(s(1:end-5), '/Engine'));
end

%  Load coil data
load coil; load coilCAD;
clear pointsXZ;

%%  Coil parameters
%  Define dIdt (for electric field)
dIdt = 9.4e7;       %   Amperes/sec (2*pi*I0/period)
%  Define I0 (for magnetic field)
I0 = 5e3;       %   Amperes

%  Parameters
mu0     = 1.25663706e-006;  %   magnetic permeability of vacuum(~air)
levels  = 20;               %   number of levels in contour plot
comp    = 2;                %   field component to be plotted (1, 2, 3 or x, y, z) 
temp        = ['x' 'y' 'z'];
label       = temp(comp);

%   Plane window (from xmin to xmax and from zmin to zmax)
scale = 1.2;
Xmin = min(strcoil.Pwire(:, 1));
Xmax = max(strcoil.Pwire(:, 1));
Zmin = min(strcoil.Pwire(:, 3));
Zmax = max(strcoil.Pwire(:, 3));
xmin = scale*Xmin;
xmax = scale*Xmax;
zmin = Zmin - 0.05; %   5 cm down
zmax = Zmax + (Zmax-Zmin);
Y    = 0.0;          %  position of the XZ plane

%   Plot the plane
f1 = figure;
bemf1_graphics_coil_CAD(P, t, 0);
patch([xmin xmin xmax xmax], [Y Y Y Y], [zmin zmax zmax zmin], 'c', 'FaceAlpha', 0.25);
view(10, 20);

%  Nodal points on the surface (MsxMs nodal points)      
Ms = 400;
%   Coronal plane
x = linspace(xmin, xmax, Ms);
z = linspace(zmin, zmax, Ms);
[X, Z]  = meshgrid(x, z);
pointsXZ(:, 1) = reshape(X, 1, Ms^2);
pointsXZ(:, 2) = Y*ones(1, Ms^2);
pointsXZ(:, 3) = reshape(Z, 1, Ms^2);    

%   Field at the surface (XZ or coronal plane, MsxMs nodal points)         
tic
Field        = bemf3_inc_field_electric(strcoil, pointsXZ, dIdt, mu0); 
fieldPlaneTime = toc  

%  Plot field on the surface (coronal plane)
f2          = figure;     
temp        = -Field(:, comp);
th1         = max(temp);
th2         = min(temp);
bemf2_graphics_vol_field(temp, th1, th2, levels, x, z);
xlabel('x, m'); ylabel('z, m');
colormap parula; colorbar;
title(strcat('E-field V/m, ', label, '-component in the coronal plane'));

%   Additionally, plot coil cross-section
[edges, TriP, TriM] = mt(t);
[Pi, ti, polymask, flag] = meshplaneintXZ(P, t, edges, TriP, TriM, Y);
if flag % intersection found                
    for n = 1:size(polymask, 1)
        i1 = polymask(n, 1);
        i2 = polymask(n, 2);
        line(Pi([i1 i2], 1), Pi([i1 i2], 3), 'Color', 'r', 'LineWidth', 2);
    end   
end

grid on; set(gcf,'Color','White');
axis equal; axis tight;
