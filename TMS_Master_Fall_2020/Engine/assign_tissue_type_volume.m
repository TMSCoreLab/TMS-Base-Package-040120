function [id] = assign_tissue_type_volume(obsPoints, n, Center, Indicator)
    %This function determines the tissue type that exists at user-specified
    %observation points based on proximity to model surfaces.
    % obsPoints (Mx3): Observation points whose tissue types are to be determined
    % n (Nx3): normal vectors of model facts
    % Center (Nx3): centroids of model facets
    % Indicator (N): tissue codes of model facets
    
    % id (M): tissue codes for observation points
    %Copyright WAW/SNM 2020
    %Preallocate
    distTriangles = zeros(size(obsPoints, 1), max(Indicator));
    pointInsideShell = -ones(size(obsPoints, 1), max(Indicator));
    %Construct matrices telling which surfaces enclose/don't enclose which observation points
    for j = 1:max(Indicator)
        eligibleCenters = Center(Indicator == j, :);
        eligibleN = n(Indicator == j, :);
        
        %For the current tissue, find the triangle closest to each observation point
        [nearestTriangles, distTriangles(:,j)] = knnsearch(eligibleCenters, obsPoints, 'K', 1);
        
        %Get vectors from observation points to nearest triangles
        r = eligibleCenters(nearestTriangles, :) - obsPoints;
        
        %If the dot product of [vector from point to triangle center] and [triangle normal] is positive, the point is inside the current shell.
        pointInsideShell(:,j) = dot(r, eligibleN(nearestTriangles, :), 2) > 0;
    end
    
    %Knowing which observation points are within which surface, find the tissue type at the observation point
    id = zeros(size(obsPoints, 1), 1);
    for j = 1:size(obsPoints, 1)
        %Catch the case where the observation point lies completely outside the model
        if(all(pointInsideShell(j,:) == 0))
            id(j) = 0;
            continue;
        end
        
        %The tissue type assigned to the point is the tissue type of the closest outside shell
        insideIndices = find(pointInsideShell(j,:));
        tempDistances = distTriangles(j, insideIndices);
        [~,insideIndicesIndex] = min(tempDistances);
        id(j) = insideIndices(insideIndicesIndex);
    end
end