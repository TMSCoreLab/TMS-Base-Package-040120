function [contrast, condin, condout] = assign_initial_conductivities(cond, condambient, Indicator, enclosingTissueIdx)
%   This function assigns conductivity contrasts to tissue boundaries.
%   cond: Interior conductivity of a given tissue
%   condambient: conductivity of the medium exterior to the outermost tissue
%   Indicator: for each facet, stores the identifier of the tissue to which that facet belongs
%   enclosingTissueIdx: for each tissue, stores the identifier of that tissue's assumed exterior tissue

%   Copyright WAW/SNM 2020

    Condout = zeros(1, length(cond));
    Condin = zeros(1, length(cond));
    Contrast = zeros(1, length(cond));
    %For each tissue, assign the default conductivity inside and outside
    for j=1:length(cond)
        if(enclosingTissueIdx(j) == 0)
            Condout(j) = condambient;
        else
            Condout(j) = cond(enclosingTissueIdx(j));
        end

        Condin(j) = cond(j);
        Contrast(j) = (Condin(j) - Condout(j))/(Condin(j) + Condout(j));
    end
    Contrast(isnan(Contrast)) = 0; %Air/Air boundaries produce NaNs

    %   Define array of contrasts for every triangular facet 
    contrast = zeros(size(Indicator, 1), 1);
    condin   = zeros(size(Indicator, 1), 1);
    condout  = zeros(size(Indicator, 1), 1);

    for m = 1:length(Contrast)
        contrast(Indicator==m)  = Contrast(m);
        condin(Indicator==m)    = Condin(m);
        condout(Indicator==m)   = Condout(m);
    end
end