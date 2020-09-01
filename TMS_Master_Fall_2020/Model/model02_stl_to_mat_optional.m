%   This script does stl to MATLAB format conversion. Also, it conditions
%   the original mesh if necessary.
%
%   Copyright SNM 2017-2020

clear all; %#ok<CLALL>
if ~isunix
    s = pwd; addpath(strcat(s(1:end-6), '\Engine'));
else
    s = pwd; addpath(strcat(s(1:end-6), '/Engine'));
end

FileName_ = uigetfile('*.stl','Select the tissue mesh file to open', 'MultiSelect', 'on');
if ~iscell(FileName_)
    FileName{1} = FileName_;
else
    FileName   = FileName_;
end
for m = 1:length(FileName)
    [P, t, normals, dummy] = stlReadAscii(FileName{m}); 
    display(FileName{m});
    NewName =  strcat(FileName{m}(1:end-4), '.mat');
    save(NewName, 'P', 't', 'normals');           
end
