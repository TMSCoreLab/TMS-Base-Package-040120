%   This script accurately computes and displays magnetic field sampled
%   over a surface (transverse plane) due to a coil via the plain FMM method
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
clear pointsXY;

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
Xmin = min(strcoil.Pwire(:, 1));
Xmax = max(strcoil.Pwire(:, 1));
Ymin = min(strcoil.Pwire(:, 2));
Ymax = max(strcoil.Pwire(:, 2));
xmin = scale*Xmin;
xmax = scale*Xmax;
ymin = scale*Ymin;
ymax = scale*Ymax ;
Z    = min(strcoil.Pwire(:, 3))-0.04;          %  position of the XY plane (4 cm down)

%   Plot the plane
f1 = figure;
bemf1_graphics_coil_CAD(P, t, 0);
patch([xmin xmin xmax xmax], [ymin ymax ymax ymin], [Z Z Z Z], 'c', 'FaceAlpha', 0.25);
view(10, 20);
        
%  Nodal points on the surface (MsxMs nodal points)      
Ms = 250;
%   Coronal plane
x = linspace(xmin, xmax, Ms);
y = linspace(ymin, ymax, Ms);
[X, Y]  = meshgrid(x, y);
pointsXY(:, 1) = reshape(X, 1, Ms^2);
pointsXY(:, 2) = reshape(Y, 1, Ms^2);    
pointsXY(:, 3) = Z*ones(1, Ms^2);  

%   Field on the surface (XY or transverse plane, MsxMs nodal points)         
tic
Field        = bemf3_inc_field_magnetic(strcoil, pointsXY, I0, mu0); 
fieldPlaneTime = toc  

%  Plot field on the surface (transverse plane)
f2      = figure;
temp    = I0*Field(:, comp);
th1     = max(temp);             %   threshold for magnetic-field plot, T
th2     = min(temp);             %   threshold for magnetic-field plot, T
bemf2_graphics_vol_field(temp, th1, th2, levels, x, y);
xlabel('x, m'); ylabel('z, m');
colormap parula; colorbar;
title(strcat('B-field T, ', label, '-component in the transverse plane'));

%   Additionally, plot coil cross-section
[edges, TriP, TriM] = mt(t);
[Pi, ti, polymask, flag] = meshplaneintXY(P, t, edges, TriP, TriM, Z);
if flag % intersection found                
    for n = 1:size(polymask, 1)
        i1 = polymask(n, 1);
        i2 = polymask(n, 2);
        line(Pi([i1 i2], 1), Pi([i1 i2], 2), 'Color', 'r', 'LineWidth', 2);
    end   
end

grid on; set(gcf,'Color','White');
axis equal; axis tight;
