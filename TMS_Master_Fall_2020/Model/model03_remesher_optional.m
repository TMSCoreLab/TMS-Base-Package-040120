%   This sript converts the stl model files to MATLAB data files and/or
%   remeshes the head model (for remeshing creates both STL and MAT data
%   files) given the maximum edge length while keeping the surface topology
%   as close as possible to the original one. The units of the original stl
%   file are mm. Note: remeshing takes time! Requires MATLAB 2019 or newer
%
%   Copyright SNM 2019-2020

clear all; %#ok<CLALL>
if ~isunix
    s = pwd; addpath(strcat(s(1:end-6), '\Engine'));
else
    s = pwd; addpath(strcat(s(1:end-6), '/Engine'));
end

MEL         = 1.0e-3;         %    maximum edge length in m (for remesher);
MELinfo     = '_1mm';         %    for data files 

FName = uigetfile('*.stl','Select the tissue mesh file to open', 'MultiSelect', 'on');
if iscell(FName)
    FileName = FName;
else
    FileName{1} = FName;
end

for m = 1:length(FileName)
    display(FileName{m});
    p = platform; p.FileName = FileName{m}; M = mesh(p,'MaxEdgeLength', MEL); [P, t] = exportMesh(p);
    P = P*1e3;                  %   save the data in mm
    [P, t]      = fixmesh(P, t(:, 1:3));
    TR          = triangulation(t(:, 1:3), P);
    stlwrite(TR, strcat(FileName{m}(1:end-4), MELinfo, '.stl'));   
    %   Finding vertices of triangles
    vert1 = P(t(:, 1),:);
    vert2 = P(t(:, 2),:);
    vert3 = P(t(:, 3),:);
    %   Finding edges
    edge1 = vert2 - vert1;
    edge2 = vert3 - vert1;
    %   Finding outer normal vectors
    normal = cross(edge1, edge2, 2);                                 %      Calculating the normal vector of a triangle
    Length = sqrt(normal(:,1).^2+normal(:,2).^2+normal(:,3).^2);     %      Calculating length of the normal vectror
    normals= normal./(repmat(Length,size(normal,3),3));              %      Normalization
    %   Saving the result
    NewName  =  strcat(FileName{m}(1:end-4), MELinfo, '.mat');
    save(NewName, 'P', 't', 'normals'); 
end
