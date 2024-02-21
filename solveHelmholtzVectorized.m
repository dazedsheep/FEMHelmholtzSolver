function [U] = solveHelmholtzVectorized(elements, kappa, beta, f, hI, g, n)

hVec = zeros(n,1);
FVec = zeros(n,1);

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

% we need the inverse of these matrices
C = pageinv(B);

% determinant of 3x3 matrix vectorized over all triangles
determinant = n1x .* n2y + n2x .* n3y  + n3x .* n1y - n1x .* n3y - n2x .* n1y - n3x .*n2y;
area = abs(determinant);

% calculate stiffness matrix
% the gradients of the P1 basis functions can be pre-computed
basisFuncGradients = [-1 -1; 1 0; 0 1];

% basisFuncGradients * inv(B(:,:,1)) * (basisFuncGradients * inv(B(:,:,1)))' * area(1)/2
K1 = pagemtimes(basisFuncGradients,C);

% finally, the stiffness matrix
K = pagemtimes(pagemtimes(K1,pagetranspose(K1)), reshape(area.*1/2, 1, 1, size(area,2)));

% calculate mass matrix
M_t = area/24 .* [2; 1; 1; 1; 2; 1; 1; 1; 2];

% local to global indices
rowK = elements.tri(:, [1 2 3 1 2 3 1 2 3]).';
colK = elements.tri(:, [1 1 1 2 2 2 3 3 3]).';

Z_t = reshape(K, 9, size(K,3)) - kappa^2 * M_t;
% put the elements in Z_t on the right position with the local to global
% index
Z = sparse(rowK, colK, Z_t, size(elements.points,1),size(elements.points,1));

%% process the transparent (Neumann or Robin) boundary conditions
e_Vec = elements.points(elements.bedges(:,1),:) - elements.points(elements.bedges(:,2),:);
e_len = sqrt(sum(e_Vec.^2,2));
t_bM = e_len'/6 .* [2;1;1;2]; % boundary element mass matrix
brow = elements.bedges(:,[1,2,1,2]);
bcol = elements.bedges(:,[1,1,2,2]);
tBM = sparse(brow,bcol, t_bM, size(elements.points,1),size(elements.points,1));

M = sparse(rowK,colK, M_t, size(elements.points,1),size(elements.points,1));
A = Z + 1i*beta*kappa*tBM;

% get the right hand side
%FVec(unique(elements.nodeIndex(:,1)),1) = f(unique(elements.nodeIndex(:,1)));

% get the neumann/robin boundary values
%hVec(elements.bedges(:,1)) = hI(elements.bedges(:,1));

% right hand side
b = tBM*hI + M*f;

%% solve the system
U = A\b;

end

