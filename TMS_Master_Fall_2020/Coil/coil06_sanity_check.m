%   This script performs magnetic field computations for a single volumetric
%   ring with M wires and compares the result with the analytical
%   solution
%
%   Copyright SNM 2018-2020

clear all; %#ok<CLALL>
if ~isunix
    s = pwd; addpath(strcat(s(1:end-5), '\Engine'));
else
    s = pwd; addpath(strcat(s(1:end-5), '/Engine'));
end

%   The coil includes one single conductor in the form of a ring. 
%   The conductor centerline model is given first
theta = [0:pi/50:2*pi];
r = 0.02;                          %   ring radius in m
x = r.*cos(theta);                 %   ring
y = r.*sin(theta);                 %   ring

%   Other parameters
a    = 0.002;       %   conductor diameter in m
M    = 16;          %   number of cross-section subdivisions 
flag = 1;           %   circular cross-section    
sk   = 1;           %   surface current distribution (skin layer 1, uniform 0)

%   Create CAD and wire models for the single conductor
Pcenter(:, 1) = x';
Pcenter(:, 2) = y';
Pcenter(:, 3) = 0;
strcoil       = meshwire(Pcenter, a, a, M, flag, sk);
[P, t]        = meshsurface(Pcenter, a, a, M, flag);  %   CAD mesh (optional, slow)     
tind          = ones(size(t, 1), 1); 

%   Display CAD and wire models
f1 = figure;
bemf1_graphics_coil_CAD(P, t, 0);
%bemf1_graphics_coil_wire(strcoil, [0 1 0]); 
view(-4, 24);

%%  Define nodal points along the observation line (1xM nodal points)
N = 51;        
zline      = linspace(-0.1, +0.1, N);
pointsline(1:N, 1) = 0;
pointsline(1:N, 2) = 0;
pointsline(1:N, 3) = zline';  

%%  Compute the magnetic field via the FMM given 1 A of current
mu0           = 1.25663706e-006;  %   Magnetic permeability of vacuum(~air)
Binc          = bemf3_inc_field_magnetic(strcoil, pointsline, 1, mu0); 

%%  Compute the magnetic field via the direct integration (Biot-Savart law)
segvector   = Pcenter(2:end, :) - Pcenter(1:end-1, :);
segmoment   = sqrt(dot(segvector, segvector, 2));
moments     = segmoment;
segpoints   = 0.5*(Pcenter(2:end, :) + Pcenter(1:end-1, :));   
directions  = segvector./repmat(segmoment, 1, 3);

tic
Binc1       = zeros(size(pointsline, 1), 3);      %   standard format 
%   B(m, :) = +mu0/(4*pi)*...
%   sum(moments(n)*cross[directions(n, 3), (points(m, :) - segpoints(n, :))]/norm(points(m, :) - segpoints(n, :))^3)
for m =1:size(pointsline, 1)
        temp        = repmat(pointsline(m, :), size(segpoints, 1), 1) - segpoints;      %   these are distances to the observation point
        DIST        = sqrt(dot(temp, temp, 2));                                         %   single column        
        I           = repmat(moments, 1, 3).*cross(directions, temp, 2)./repmat(DIST.^3, 1, 3);
        Binc1(m, :)     = mu0/(4*pi)*sum(I);
end
DirectTime = toc

%%  Compute the magnetic field via the analytical formula
Binc2       = mu0/2*r^2./(r^2 + zline.^2).^(3/2);

%%  Plot all fields
f2 = figure;
plot(zline, Binc(:, 3), '-r', zline, Binc1(:, 3), '*b', zline, Binc2, 'om');
grid on; 
xlabel('distance along the loop axis, m');
ylabel('Magnetic field strength in T given 1 A of current');
title('Comparison with the analytical solution for a single loop: on axis B field in T')




