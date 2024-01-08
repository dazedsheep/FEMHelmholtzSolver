function [U] = solveHelmholtzCVectorized(elements, omega, kappa, gamma, beta, f, hI, g, n)
% we fix the quadrature so we can unroll the integration
%N = 2;
%[quadratureParameters.Points, quadratureParameters.W] = triangleQuadrature(N);
quadratureParameters.Points = [0.280019915499074,0.644948974278318;0.666390246014701,0.155051025721682;0.075031110222608,0.644948974278318;0.178558728263616,0.155051025721682].';
quadratureParameters.W = [0.090979309128011;0.159020690871989;0.090979309128011;0.159020690871989];

hVec = zeros(n,1);

phi1 = @(x,y) 1-x-y;
phi2 = @(x,y) x;
phi3 = @(x,y) y;

bfunc1 = @(x,y) phi1(x,y).*phi1(x,y);
bfunc2 = @(x,y) phi1(x,y).*phi2(x,y);
bfunc3 = @(x,y) phi1(x,y).*phi3(x,y);
bfunc4 = @(x,y) phi2(x,y).*phi2(x,y);
bfunc5 = @(x,y) phi2(x,y).*phi3(x,y);
bfunc6 = @(x,y) phi3(x,y).*phi3(x,y);

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
% 
% we need to transform the quadrature points for all the elements
sPo = repmat(reshape([n1x;n1y], 2,1,size(elements.triangles,2)), 1,4,1);
quP = repmat(quadratureParameters.Points,1,1,size(elements.triangles,2));
tQuP = pagemtimes(B,quP) + sPo;

MC = zeros(3,3,size(elements.triangles,2));

KappaSq = double(kappa(tQuP(1,:,:),tQuP(2,:,:)).^2);
A1 = bfunc1(quP(1,:,:), quP(2,:,:));
A2 = bfunc2(quP(1,:,:), quP(2,:,:));
A3 = bfunc3(quP(1,:,:), quP(2,:,:));
A4 = bfunc4(quP(1,:,:), quP(2,:,:));
A5 = bfunc5(quP(1,:,:), quP(2,:,:));
A6 = bfunc6(quP(1,:,:), quP(2,:,:));
weights = repmat(quadratureParameters.W.',1,1,size(elements.triangles,2));
MC(1,1,:) = squeeze(pagemtimes(weights, pagetranspose(A1(:,:,:).*KappaSq(:,:,:)))).'.*area;
MC(1,2,:) = squeeze(pagemtimes(weights, pagetranspose(A2(:,:,:).*KappaSq(:,:,:)))).'.*area;
MC(1,3,:) = squeeze(pagemtimes(weights, pagetranspose(A3(:,:,:).*KappaSq(:,:,:)))).'.*area;

MC(2,1,:) = MC(1,2,:);
MC(2,2,:) = squeeze(pagemtimes(weights, pagetranspose(A4(:,:,:).*KappaSq(:,:,:)))).'.*area;
MC(2,3,:) = squeeze(pagemtimes(weights, pagetranspose(A5(:,:,:).*KappaSq(:,:,:)))).'.*area;

MC(3,1,:) = MC(1,3,:);
MC(3,2,:) = MC(2,3,:);
MC(3,3,:) = squeeze(pagemtimes(weights, pagetranspose(A6(:,:,:).*KappaSq(:,:,:)))).'.*area;

% local to global indices
rowK = elements.tri(:, [1 2 3 1 2 3 1 2 3]).';
colK = elements.tri(:, [1 1 1 2 2 2 3 3 3]).';

Z_t = reshape(K, 9, size(K,3))  - reshape(MC,9,size(MC,3));
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
U = gather(A\b);

end

