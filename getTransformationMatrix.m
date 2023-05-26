function [A] = getTransformationMatrix(triangleVertices)
% transformation matrix from the unit triangle to the actual element
A = [triangleVertices(1,2) - triangleVertices(1,1), triangleVertices(1,3) - triangleVertices(1,1);triangleVertices(2,2) - triangleVertices(2,1), triangleVertices(2,3) - triangleVertices(2,1)];
end

