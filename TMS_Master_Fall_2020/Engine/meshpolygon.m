%   This function assembles a polygon from a list of arbitrary edges
%
%   Copyright SNM 2012-2020
function [Ppolygon] = meshpolygon(P, E)
    E      = sortrows(E);
    stop   = size(E, 1); 
    nextedge = 1;
    node = E(1, 2);
    for m = 1:stop
        ind1 = find(E(:, 1)==node);
        ind2 = find(E(:, 2)==node);
        ind1 = setdiff(ind1, nextedge);
        ind2 = setdiff(ind2, nextedge);
        if ind1 nextedge = [nextedge ind1]; node = E(ind1, 2); end
        if ind2 nextedge = [nextedge ind2]; node = E(ind2, 1); end
    end
    E = E(nextedge, :);
    PS   = E(1, :);
    stop = size(E, 1);
    for m = 2:stop
        PS = [PS setdiff(E(m, :), PS)];
    end
    PS = [PS 1];
    Ppolygon = P(PS, :);
end