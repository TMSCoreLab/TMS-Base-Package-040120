function P = meshrotate2(P, axis, theta)
%   SYNTAX
%   P = meshrotate2(P, axis, theta)
%   DESCRIPTION
%   This function implements rotation (Rodrigues' rotation formula) of the
%   mesh about an axis with vector 'axis'. Rotation angle is to be given in
%   radians.
%
%   To display the mesh use: fv.faces = t; fv.vertices = P;
%   patch(fv, 'FaceColor', 'y'); axis equal; view(160, 60); grid on;
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

    NX   = axis(1)/norm(axis);        %   unit vector
    NY   = axis(2)/norm(axis);        %   unit vector
    NZ   = axis(3)/norm(axis);        %   unit vector
    k    = [NX NY NZ]';
    K     = repmat(k, [1 size(P, 1)])'/sqrt(dot(k, k));
    P = P*cos(theta) + cross(K, P, 2)*sin(theta) + K.*repmat(dot(K, P, 2), [1 3])*(1-cos(theta));
end