%   "Coil tester line electric" tool for the entire coil. Outputs
%   solenoidal electric field -dA/dt (or any of the field components) along
%   the line given I0 A of conductor current and frequency freq

%   Copyright SNM 2017-2020

clear all; %#ok<CLALL>
if ~isunix
    s = pwd; addpath(strcat(s(1:end-5), '\Engine'));
else
    s = pwd; addpath(strcat(s(1:end-5), '/Engine'));
end
load coil; load coilCAD;

%%  Define EM constants
mu0         = 1.25663706e-006;  %   Magnetic permeability of vacuum(~air)

%%  Coil parameters
%  Define dIdt (for electric field)
dIdt = 9.4e7;       %   Amperes/sec (2*pi*I0/period)
%  Define I0 (for magnetic field)
I0 = 5e3;       %   Amperes

%%  Define nodal points along the line (1xM nodal points)
M = 1001;        
line      = linspace(-0.05, 0.05, M);
pointsline(1:M, 1) = 0;
pointsline(1:M, 2) = 0;
pointsline(1:M, 3) = line';   

%%   Plot the line
f1 = figure;
hold on;
bemf1_graphics_coil_CAD(P, t, 0);
plot3(pointsline(:, 1), pointsline(:, 2), pointsline(:, 3), '-r', 'lineWidth', 2);
view(10, 20);

%%  Find the E-field on the line        
tic
Einc           = bemf3_inc_field_electric(strcoil, pointsline, dIdt, mu0); 
fieldTime = toc  

%%  Graphics for the line
f2 = figure;
hold on;
plot(line, +Einc(:, 1), '-r', 'LineWidth', 2); 
plot(line, +Einc(:, 2), '-m', 'LineWidth', 2); 
plot(line, +Einc(:, 3), '-b', 'LineWidth', 2); 
grid on;

title('Line field E in Vm, red:Ex magenta:Ey blue:Ez');
xlabel('Distance from the origin, m');
ylabel('Electric field, V/m');
set(gcf,'Color','White');