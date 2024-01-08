function [U] = solveHelmholtzCVectorizedKappaSampled(elements, omega, kappa, gamma, beta, f, hI, g, n)

hVec = zeros(n,1);

% initialize the needed data
% extract nodes numbers of the 3 vertices of each triangle 
n1x = elements.points(elements.tri(:,1),1).'; 
n1y = elements.points(elements.tri(:,1),2).'; 
n2x = elements.points(elements.tri(:,2),1).'; 
n2y = elements.points(elements.tri(:,2),2).'; 
n3x = elements.points(elements.tri(:,3),1).'; 
n3y = elements.points(elements.tri(:,3),2).';  
m = size(elements.tri,1);
n1 = [n1x;n1y];
n2 = [n2x;n2y];
n3 = [n3x;n3y];

% compute the element wise transformation matrices
A = [n2(1,:) - n1(1,:), n3(1,:) - n1(1,:); n2(2,:) - n1(2,:), n3(2,:) - n1(2,:)];
B(1,1,:) = A(1,1:m);
B(1,2,:) = A(1,(m+1):2*m);
B(2,1,:) = A(2,1:m);
B(2,2,:) = A(2,(m+1):2*m);

% we need the inverse of each these matrices
C = pageinv(B);

% determinant of 3x3 matrix vectorized over all triangles
determinant = n1x .* n2y + n2x .* n3y  + n3x .* n1y - n1x .* n3y - n2x .* n1y - n3x .*n2y;
area = abs(determinant);

% calculate stiffness matrix
% the gradients of the P1 basis functions can be pre-computed
basisFuncGradients = [-1 -1; 1 0; 0 1];

% basisFuncGradients * inv(B(:,:,1)) * (basisFuncGradients * inv(B(:,:,1)))' * area * 1/2
K1 = pagemtimes(basisFuncGradients,C);

% finally, the stiffness matrix
K = pagemtimes(pagemtimes(K1,pagetranspose(K1)), reshape(area.*1/2, 1, 1, size(area,2)));

% calculate mass matrix without spatial coefficient
M_t = area/24 .* [2; 1; 1; 1; 2; 1; 1; 1; 2];

% calculate weighted mass matrix
% this just takes one value for kappa for each of the triangles and
% quadrature points

KappaSq = repmat(mean(kappa(elements.nodeIndex).^2,2),1,9,1).';
% kappaPointsIdx = elements.nodeIndex(:,1);
% KappaSq = repmat(kappa(kappaPointsIdx).^2,1,9,1).';
MC = M_t .* KappaSq;

% local to global indices
rowK = elements.tri(:, [1 2 3 1 2 3 1 2 3]).';
colK = elements.tri(:, [1 1 1 2 2 2 3 3 3]).';

Z_t = reshape(K, 9, size(K,3))  - MC;
% put the elements in Z_t on the right position with the local to global
% index
Z = sparse(rowK, colK, Z_t, size(elements.points,1),size(elements.points,1));

% compute length of boundary edges
e_Vec = elements.points(elements.bedges(:,1),:) - elements.points(elements.bedges(:,2),:);
e_len = sqrt(sum(e_Vec.^2,2));

% boundary element mass matrix
t_bM = e_len'/6 .* [1;2;2;1]; 

% local to global index for the boundary
brow = elements.bedges(:,[1,2,1,2]);
bcol = elements.bedges(:,[1,1,2,2]);

% sparse boundary mass matrix
tBM = sparse(brow,bcol, t_bM, size(elements.points,1),size(elements.points,1));

% sparse mass matrix
M = sparse(rowK,colK, M_t, size(elements.points,1),size(elements.points,1));

A = Z + (1i*beta*omega+gamma)*tBM;

% get the neumann/robin boundary values
hVec(elements.bedges(:,1)) = hI(elements.bedges(:,1));

% right hand side
b = tBM*hVec - M*f;

%% solve the system
U = A\b;

end