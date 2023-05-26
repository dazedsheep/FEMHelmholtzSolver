function [val] = numIntTriangle(funcIn, TriangleVertices, quadratureParameters)
% This function uses the quadrature rule on the unit triangle for numerical integration
% 
% f                     ... function to be integrated
% TriangleVertices      ... vertices of triangle
% quadratureParameters  ... weights and evaluation points for gaussian
% quadrature on the unit triangle

% the baricentric points are transformed from the unit triangle to the
% actual reference triangle element
transformMatrix = getTransformationMatrix(TriangleVertices);

[X,Y] = transformFromUnitTriangle(quadratureParameters.Points, TriangleVertices);

fvals = zeros(size(X,1),1);

for i=1:size(X,1)
    fvals(i) = feval(funcIn,X(i),Y(i));
end

val = quadratureParameters.W' * fvals * det(transformMatrix);

end