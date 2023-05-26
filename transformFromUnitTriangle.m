function [tX,tY] = transformFromUnitTriangle(Points, TriangleVertices)
% X and Y are coordinates in the unit triangle which shall be converted to
% the triangle given by TriangleVertices.

% X and Y have the same number of elements
tX = zeros(size(Points,1),1);
tY = zeros(size(Points,1),1);
tM = getTransformationMatrix(TriangleVertices);
for i=1:size(Points,1)
%   for j=1:size(X,2)
        v = tM*[Points(i,1);Points(i,2)] + TriangleVertices(:,1);
        tX(i) = v(1);
        tY(i) = v(2);
%    end
end
end

