function [K, M, F] = assembleFEM(element, f, quadratureParameters)

basisFuncDiff = [-1 -1; 1 0; 0 1];
basisFunctions = 3;

F = zeros(3,1);

% the transformation matrix for the specific element
tM = getTransformationMatrix(element.triangleVertices);
% its inverse
tMI = inv(tM);
% and its determinant
dtM= det(tM);

K = zeros(basisFunctions, basisFunctions);
F = zeros(basisFunctions,1);
% build the stiffness element matrix 
for i=1:basisFunctions

    for j=1:basisFunctions  
            
            dotProductNablaPhis = basisFuncDiff(i,:)*tMI * (basisFuncDiff(j,:)*tMI)';
            % we use triangle quadrature to compute the element-wise
            % integrals
            K(i,j) = numIntTriangle(@(x,y) dotProductNablaPhis, element.triangleVertices, quadratureParameters);

    end

end
    
% the mass element matrix does not change as the basis functions are
% fixed, it is just scaled by the determinant of the transformation
% matrix
M = dtM .* 1/24 .* [2 1 1 ; 1 2 1 ; 1 1 2];

% get the right hand sinde
F(1) = f(element.triangleVertices(1,1),element.triangleVertices(2,1));
F(2) = f(element.triangleVertices(1,2),element.triangleVertices(2,2));
F(3) = f(element.triangleVertices(1,3),element.triangleVertices(2,3));

end
