%   This script computes and displays tissue intersection points (in
%   meters) with an arbitrary ray/segment
%
%   Copyright SNM 2017-2020

clear all; %#ok<CLALL>
if ~isunix
    s = pwd; addpath(strcat(s(1:end-6), '\Engine'));
else
    s = pwd; addpath(strcat(s(1:end-6), '/Engine'));
end

%   Ray parameters
orig = [0 0 0];     %   ray origin
dir  = [0 0 1];     %   ray direction
dist = 1000;        %   ray length (finite segment, in mm here)

FileName = uigetfile('*.mat','Select the tissue mesh file to open'); load(FileName, '-mat');

S = load(FileName);
d = meshsegtrintersection(orig, dir, dist, S.P, S.t);
IntersectionPoint = d(find(d>0))

