function [] = bemf2_graphics_base(P, t, c)
%   Surface plot

%   Copyright SNM 2017-2020

    p = patch('vertices', P, 'faces', t);
    p.FaceColor = c.FaceColor;
    p.EdgeColor = c.EdgeColor;
    p.FaceAlpha = c.FaceAlpha;
    daspect([1 1 1]);      
    xlabel('x, m'); ylabel('y, m'); zlabel('z, m');
	
    NumberOfTrianglesInShell = size(t, 1);
    edges = meshconnee(t);
    temp = P(edges(:, 1), :) - P(edges(:, 2), :);
    AvgEdgeLengthInShell = mean(sqrt(dot(temp, temp, 2)));
end