function d = meshsegtrintersection(orig0, dir0, dist0, P0, t0)
%   This function checks whether or not a segment characterized by orig0,
%   dir0, dist0 intersects a manifold mesh P0, t0.
%   Inputs:
%   orig0   - Origin of the segment (1 x 3)
%   dir0    - Normalized direction of the segment from origin (1 x 3)
%   dist0   - Length of the segment (1 x 1)
%   P0, t0  - Triangulation to be tested 
%   Output:
%   d -  Distances of points of intersection from the origin of the segment
%   (N x 1). If there is no intersection, the corresponding field is zero.
%   The tolerance is given internally
%   The script implements the method described in 
%   Tomas Moeller and Ben Trumbore, “Fast, Minimum Storage Ray/Triangle
%   Intersection”, Journal of Graphics Tools, 2(1):21—28, 1997
%   See also
%   http://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
%   http://www.mathworks.com/matlabcentral/fileexchange/33073-triangle-ray-intersection
%   Authors: Vishal Rathi (vkrathi@wpi.edu)
%   Janakinadh Yanamadala (jyanamadala@wpi.edu), SNM (makarov@wpi.edu)
%
%   To display the mesh use: fv.faces = d; fv.vertices = P;
%   patch(fv, 'FaceColor', 'y'); axis equal; view(160, 60); grid on;

    vert1 = P0(t0(:, 1),:);
    vert2 = P0(t0(:, 2),:);
    vert3 = P0(t0(:, 3),:);
    orig = repmat(orig0, size(vert1, 1),1);             
    dist = repmat(dist0, size(vert1, 1),1);
    dir  = repmat(dir0, size(vert1, 1),1);

    % Initialization of u,v and d
    u = zeros (size(vert1,1),1);
    d = u; v = u;

    % Finding edges
    edge1 = vert2 - vert1;
    edge2 = vert3 - vert1;

    tvec = orig - vert1;                            %   Distance to vert1 from segment origin
    pvec = cross(dir, edge2, 2);                    %   Parameter to calculate u
    det  = dot(edge1, pvec, 2);                     %   Determinant of matrix M
    parallel = abs(det)< 1024*eps*max(abs(det));    %   To test edges parallel with the segment
    if all(parallel)                                %   If all parallel then no intersections
        return;
    end

    det(parallel) = 1;              %   To avoid division by zero
    inv_det = 1.0 ./ det;           %   Find inverse of the determinant
    u = dot(tvec,pvec,2);           %   Calculate the u parameter
    u = u.*inv_det;

    % Conditional tests for u and v
    layer1 = (~ parallel & u<0 | u>1);
    if all(layer1)
        return;
    end

    qvec (~layer1,:) = cross(tvec(~layer1,:), edge1(~layer1,:), 2);             %   Parameter to calculate v
    v (~layer1,:) = dot(dir(~layer1,:),qvec(~layer1,:),2).*inv_det(~layer1,:);  %   Calculate v
    layer2 = (v<=0 | u+v>1);
    if all(layer2)
        return;
    end

    layer = (~layer1&~layer2);
    d(layer,:) = dot(edge2(layer,:),qvec(layer,:),2).*inv_det(layer,:);         %   Calculate d
    d(d<0 | d>dist) = 0;                                                        %   Compare distances and d
    d(parallel) = 0;                                                            %   Avoid values of d in parallel cases
    d(isnan(d))= 0;                                                             %   Avoid NaN (Not-a-Number) when the right-angled triangles are present
end


