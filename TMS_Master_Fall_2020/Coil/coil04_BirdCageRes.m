%   "Wire coil creator" tool for the entire birdcage coil. Models the outer surface
%   of coil windings by multiple wire segments. The output is saved in the
%   binary file coil.mat and includes:
%   Pwire{:, 3} - set of nodes for all wires 
%   Ewire{:, 2} - set of edges for all wires (default current flows from the first 
%   to the second node)
%   Swire{:, 1} - current strength weight for every segment/edge (from -1 to 1)
%   This parameter depends on the number of wires M in the indivudual
%   conductor and scales by default as 1/M.  

%   Copyright SNM 2020

clear all; %#ok<CLALL>
s = pwd; addpath(strcat(s(1:end-5), '\Engine'));

% Coil geometry definitions
IN          = 0.0254;       % inch to meter conversion factor
r_coil      = 18.52*IN;     % Radius from coil centerline to endring/rung centerline
l_rung      = 45.25*IN;     % Rung length
d_rung      = 0.375*IN;     % Rung tube outer diameter
d_endring   = 1.375*IN;     % Endring tube outer diameter
n_rungs     = 66; % 144     % Number of rungs

%  Define wire cross-section
type = 'circle';
a   = d_endring/2;         %   radius (for circle) or x-side, m  (for a rectangle/tape cross-section)
b   = d_endring/2;         %   radius (for circle) or y-side, m  (for a rectangle/tape cross-section)
M   = 12;                  %   number of cross-section subdivisions (wires)
%  Define the number of angular subdivisions per loop
N   = 256; 
%  Define positions of every loop within the coil (the initial coil axis is the z-axis)
%   Give the x-coordinates for the centerline of every loop in the xz-plane:
x   = r_coil * [1 1];   %   centerline intersection with the xz-plane
%   Give the y-coordinates for the centerline of every loop in the yz plane
y   = r_coil * [1 1];   %   centerline intersection with the xy-plane
%   Give the z-coordinates for the centerline of every loop (the z offset)
z  = (l_rung/2+d_endring/2) * [-1 1];
K  = length(z);             %    number of loops in the coil    

%   Create the base model with K loops oriented along the z-axis by cloning one loop with M wires
PwireL          = [];           % : by 3 point array
EwireL          = [];           % : by 2 edge array
SwireL          = [];           % : by 1 indicator (current strength) array for each segment/edge
NFwireL         = ones(K*M, 1); % : by 1 array of the number of nodes in every wire

for m = 1:K                 % loop over coil loops   
    [Pwire, Ewire]  = wire_single_loop(x(m), y(m), a, b, N, M, type);  %  this is the major function 
    Pwire(:, 3)     = Pwire(:, 3) + z(m); 
    EwireL          = [EwireL; Ewire+size(PwireL, 1)];
    PwireL          = [PwireL; Pwire];
    node_x          = Pwire(:,1);
    wire_x          = mean(node_x(Ewire),2);
    node_y          = Pwire(:,2);
    wire_y          = mean(node_y(Ewire),2);
    [theta, ~]      = cart2pol(wire_x, wire_y);
    ringcurrent     = sin(theta) + j*sin(theta + pi/2); %   variable ring current
    SwireL          = [SwireL; 1/M*sign(z(m))*ringcurrent];    
    for n = 1:M
       NFwireL(M*(m-1)+n) = N;  
    end
end

%%   Create the rung structure 
%  Define wire cross-section
type = 'circle';
a   = d_rung/2;            %   radius (for circle) or x-side, m  (for a rectangle/tape cross-section)
b   = d_rung/2;            %   radius (for circle) or y-side, m  (for a rectangle/tape cross-section)
M   = 12;                  %   number of cross-section subdivisions (wires)
%  Define the number of length subdivisions along the z-axis
N   = 256; 
% Define position of every rung and rung current strength/weights
L = l_rung;                     %   length of one rung
Rungs = n_rungs;                %   number of rungs in the coil  
phi = 2*pi*(0:Rungs-1)/Rungs;   %   modified by Gene
%   Give the x-coordinates for the centerline of every rung in the xy-plane:
x   = r_coil*cos(phi);     %   centerline intersection with the xy-plane
%   Give the y-coordinates for the centerline of every rung in the xy plane
y   = r_coil*sin(phi);     %   centerline intersection with the xy-plane

%   Create the base model with "Rungs" rungs oriented along the z-axis by cloning one rung with M wires
PwireR          = [];               % : by 3 point array
EwireR          = [];               % : by 2 edge array
SwireR          = [];               % : by 1 indicator (current strength) array for each segment/edge
NFwireR         = ones(Rungs*M, 1); % : by 1 array of the number of nodes in every wire
S               = N-1;              % number of segments in a single wire 
%   Define and normalize rung current
NRungsHalf      = Rungs/2;          % modified by Gene
arg             = (0:NRungsHalf-1)/NRungsHalf;
I0              = 2/sum(sin(pi*arg));
rungcurrent     = I0*cos(phi) + j*I0*cos(phi+pi/2);

for m = 1:Rungs                 % loop over coil loops   
    [Pwire, Ewire]  = wire_single_rung(L, a, b, N, M, 'circle');  %  this is the major function     
    Pwire(:, 1)     = Pwire(:, 1) + x(m); 
    Pwire(:, 2)     = Pwire(:, 2) + y(m); 
    EwireR          = [EwireR; Ewire+size(PwireR, 1)];
    PwireR          = [PwireR; Pwire];
    SwireR          = [SwireR; 1/M*repmat(rungcurrent(m), M*S, 1)];
    for n = 1:M
       NFwireR(M*(m-1)+n) = N;  
    end
end

%%   Combine two structures together
Pwire           = [PwireL; PwireR];
Ewire           = [EwireL; EwireR+size(PwireL, 1)];
Swire           = [SwireL; SwireR];
NFwire          = [NFwireL; NFwireR];
segments        = size(Swire, 1); 

%%  Plot the coil including current strength
SwirePlot   = real(Swire);    %   one mode here
Vertices    = size(Pwire, 1);
map         = jet(Vertices);
color       = zeros(Vertices, 3);
MAX         = max(SwirePlot);
MIN         = min(SwirePlot);
for m = 1:segments  %   Simple fix but it works 
    arg = (SwirePlot(m)-MIN)/(MAX-MIN);
    arg = mypalette(arg);    
    arg = ceil(arg*segments+eps);    
    color(Ewire(m, 1), :) = map(arg, :);
    color(Ewire(m, 2), :) = map(arg, :);
end
patch('faces', Ewire, 'vertices', Pwire, 'FaceVertexCData', color, 'EdgeColor','interp', 'FaceColor','none','LineWidth', 2);
axis 'equal';  axis 'tight';      
xlabel('x, m'); ylabel('y, m'); zlabel('z, m');
view(-15, 30);  set(gcf,'Color','White');
title('Coil model with the relative current strength (-1 to 1)')
%   Insert colorbar
divisions = 20;
for m = 1:divisions+1
     mytickmap{m} = num2str(-1+(m-1)/divisions*2);
end
ticks      = linspace(0, 1, divisions+1);
ticks      = mypalette(ticks);
colormap(jet) %// apply colormap
colorbar('Ticks', ticks, 'TickLabels', mytickmap, 'Location', 'east');

%%   Save coil data
strcoil.Pwire   = Pwire;
strcoil.Ewire   = Ewire;
strcoil.Swire   = Swire;
strcoil.NFwire  = NFwire;
save('coil', 'strcoil');

% Palette with the high-resolution close to zero
function y = mypalette(x)    
    y = 0.5*sign(x-0.5).*(abs(x-0.5)/0.5).^(1/2)+0.5;
end
 
 