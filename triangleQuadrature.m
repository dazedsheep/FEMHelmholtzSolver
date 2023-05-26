function [Points,Weights] = triangleQuadrature(N)
% Gaussian-Legendre quadrature for triangles
% we calculate the quadrature points on a (unit) square and then collapse
% it to a (unit) triangle
dim = 2;

for i = 1:dim
    
    [qp{i}, w{i}] = squareQuadrature(N,dim-i);

end
[QP{1:dim}] = ndgrid(qp{:});
[W{1:dim}] = ndgrid(w{:});
q = reshape(cat(dim, QP{:}), N^dim,dim);
w = reshape(cat(dim, W{:}), N^dim,dim);
% create triangle vertices
m=dim+1; 
I = eye(m);
vert=eye(m,dim);
I(2:m, 1) = -1;
c = I*vert;
% transform points to triangle
Weights = abs(det((c(2:m,:))))*prod(w,2);
qcumprods = cumprod(q,2);
Points =  [ones(N^dim, 1) [(1-q(:,1:dim-1)), ones(N^dim,1)].*[ones(N^dim,1), qcumprods(:,1:dim-2), qcumprods(:,dim)]]*c;
end

% legendre polynomials on the square
function [x,w] = squareQuadrature(N,dim)
c1 = dim + 1;
c2 = dim + 2;
n = 1:N;    % all the orders
evn = 2*n + dim;
A = [dim/c2 repmat(dim^2,1,N)./(evn.*(evn+2))];
n = 2:N;
evn = evn(n);
nek = n + dim;
evn2 = evn.^2;
U = 4*c1/(c2^2*(dim+3));
B = 4*(n.*nek).^2./(evn2.^2-evn2);
L = [A' [(2^c1)/c1; U; B']];
h = sqrt(L(2:N,2));
[V,X] = eig(diag(L(1:N,1),0) + diag(h, -1) + diag(h,1));
[X,idx] = sort(diag(X));
x = (X+1)/2;
w = (1/2)^c1 * L(1,2) * V(1,idx)'.^2;
end