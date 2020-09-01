%   This script accurately computes and displays magnetic field sampled
%   over a surface (sagittal plane) due to a coil via the plain FMM method
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
clear pointsYZ;

%%  Coil parameters
%  Define dIdt (for electric field)
dIdt = 9.4e7;       %   Amperes/sec (2*pi*I0/period)
%  Define I0 (for magnetic field)
I0 = 5e3;       %   Amperes

%  Parameters
mu0     = 1.25663706e-006;  %   magnetic permeability of vacuum(~air)
levels  = 20;               %   number of levels in contour plot
comp    = 1;                %   field component to be plotted (1, 2, 3 or x, y, z) 
temp        = ['x' 'y' 'z'];
label       = temp(comp);

%   Plane window (from xmin to xmax and from zmin to zmax)
scale = 1.2;
Ymin = min(strcoil.Pwire(:, 2));
Ymax = max(strcoil.Pwire(:, 2));
Zmin = min(strcoil.Pwire(:, 3));
Zmax = max(strcoil.Pwire(:, 3));
ymin = scale*Ymin;
ymax = scale*Ymax;
zmin = Zmin - 0.05; %   5 cm down
zmax = Zmax + (Zmax-Zmin);
X    = 0.0;          %  position of the YZ plane

%   Plot the plane
f1 = figure;
bemf1_graphics_coil_CAD(P, t, 0);
patch([X X X X], [ymin ymin ymax ymax], [zmin zmax zmax zmin], 'c', 'FaceAlpha', 0.25);
view(10, 20);
        
%  Nodal points on the surface (MsxMs nodal points)      
Ms = 250;
%   Sagittal plane
y = linspace(ymin, ymax, Ms);
z = linspace(zmin, zmax, Ms);
[Y, Z]  = meshgrid(y, z);
pointsYZ(:, 1) = X*ones(1, Ms^2);
pointsYZ(:, 2) = reshape(Y, 1, Ms^2);
pointsYZ(:, 3) = reshape(Z, 1, Ms^2);        

%   Field on the surface (YZ or sagittal plane, MsxMs nodal points)         
tic
Field        = bemf3_inc_field_magnetic(strcoil, pointsYZ, I0, mu0); 
fieldPlaneTime = toc  

%  Plot field on the surface (sagittal plane)
f2          = figure;
temp        = I0*Field(:, comp);
th1         = max(temp);
th2         = min(temp);
bemf2_graphics_vol_field(temp, th1, th2, levels, y, z);
xlabel('y, m'); ylabel('z, m');
colormap parula; colorbar;
title(strcat('B-field T, ', label, '-component in the sagittal plane'));

%   Additionally, plot coil cross-section
[edges, TriP, TriM] = mt(t);
[Pi, ti, polymask, flag] = meshplaneintYZ(P, t, edges, TriP, TriM, X);
if flag % intersection found                
    for n = 1:size(polymask, 1)
        i1 = polymask(n, 1);
        i2 = polymask(n, 2);
        line(Pi([i1 i2], 2), Pi([i1 i2], 3), 'Color', 'r', 'LineWidth', 2);
    end   
end

grid on; set(gcf,'Color','White');
axis equal; axis tight;