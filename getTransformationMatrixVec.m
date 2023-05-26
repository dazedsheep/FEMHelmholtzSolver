function [A] = getTransformationMatrixVec(n1,n2,n3)
% transformation matrix from the unit triangle to the actual element
%A = [triangleVertices(1,2) - triangleVertices(1,1), triangleVertices(1,3) - triangleVertices(1,1);triangleVertices(2,2) - triangleVertices(2,1), triangleVertices(2,3) - triangleVertices(2,1)];

A = [n2(1,:) - n1(1,:), n3(1,:) - n1(1,:); n2(2,:) - n1(2,:), n3(2,:) - n1(2,:) ];

end

