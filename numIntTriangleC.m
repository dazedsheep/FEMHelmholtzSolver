function [val] = numIntTriangleC(a, b, X, Y, dtM, quadratureParameters)
% This function uses the quadrature rule on the unit triangle for numerical integration
% 
% f                     ... function to be integrated
% TriangleVertices      ... vertices of triangle
% quadratureParameters  ... weights and evaluation points for gaussian
% quadrature on the unit triangle

% the baricentric points are transformed from the unit triangle to the
% actual reference triangle element

%fvals = zeros(size(X,1),1);

% for i=1:size(X,1)
%     fvals(i) = feval(a,X(i),Y(i)).*feval(b, quadratureParameters.Points(i,1), quadratureParameters.Points(i,2));
% end

val = quadratureParameters.W' * double(feval(a,X,Y).*feval(b, quadratureParameters.Points(:,1), quadratureParameters.Points(:,2))) * dtM;

end