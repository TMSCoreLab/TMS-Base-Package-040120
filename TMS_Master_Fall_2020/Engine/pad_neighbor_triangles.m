function tneighbor = pad_neighbor_triangles(tneighbor)
%   If there are holes in the mesh, not all triangles will always have three neighbors.
%   In this case, duplicate the existing neighbors to fill the gaps

%   Copyright WAW/SNM 2020
    nanrows = find(isnan(tneighbor(:,1) + tneighbor(:,2) + tneighbor(:,3))); %Any operation with a NaN returns a NaN
    for j=1:length(nanrows)
        tneighbor(nanrows(j),:) = sort(tneighbor(nanrows(j),:)); %If there are any NaNs, this operation will push them to the end of their respective rows
        if(isnan(tneighbor(nanrows(j), 1)))
            tneighbor(nanrows(j),1) = nanrows(j); %If the triangle has no neighbors, make it its own neighbor
        end
        for k=2:length(tneighbor(nanrows(j),:))
            if(isnan(tneighbor(nanrows(j), k)))
                tneighbor(nanrows(j), k) = tneighbor(nanrows(j), k-1);
            end
        end
    end
    if(any(isnan(tneighbor)))
        error('tneighbor is not fully populated');
    end

end