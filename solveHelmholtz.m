function [U] = solveHelmholtz_test(elements, kappa, beta, f, hI, g, n)
%% assemble the finite elements
% todo parallelize this for loop, the elements can be easily computed in
% parallel
N = 2;
[quadratureParameters.Points, quadratureParameters.W] = triangleQuadrature(N);

idx  = 1;
basisFunctions = 3;

% initialise all the required arrays with their required size (this is
% needed for mexing code)
row  = zeros(basisFunctions^2 * size(elements.triangles,2), 1);
col  = zeros(basisFunctions^2 * size(elements.triangles,2), 1);
Zsparse = zeros(basisFunctions^2 * size(elements.triangles,2), 1);
Msparse = zeros(basisFunctions^2 * size(elements.triangles,2), 1);

hVec = zeros(n,1);
FVec = zeros(n,1);
j = 1;

for i=1:size(elements.triangles,2)
    % prepare the prescribed h on the boundary
    [K, M, F] = assembleFEM(elements.triangles{i}, f, quadratureParameters);
    
    % local index
    idx = (i-1)*basisFunctions^2 + 1;
    localNodeIdx        = idx:(idx+basisFunctions^2-1);             
    
    % global index
    globalNodeIdx       = elements.nodeIndex(i,:);  

    % local to global index (rows and columns) needed later on for the
    % sparse matrices
    row(localNodeIdx)   = globalNodeIdx([1;1;1],:); 
    glblNIdxt = globalNodeIdx';
    col(localNodeIdx)   = glblNIdxt(:,[1 1 1]);
    
    FVec(globalNodeIdx) = F;
 
    Z = K - kappa^2 * M;
    Zsparse(localNodeIdx) = Z;
    Msparse(localNodeIdx) = M;
    
end

% get the right hand side of the neumann/robin boundary conditions
for i=1:size(elements.bedges, 1)
    hVec(elements.bedges(i,1)) = hI(elements.points(elements.bedges(i,1), 1), elements.points(elements.bedges(i,1), 2));
end


%% process the transparent (Neumann or Robin) boundary conditions
e_Vec = elements.points(elements.bedges(:,1),:) - elements.points(elements.bedges(:,2),:);
e_len = sqrt(sum(e_Vec.^2,2));
t_bM = e_len'/6 .* [1;2;2;1]; % boundary element mass matrix
brow = elements.bedges(:,[1,2,1,2]);
bcol = elements.bedges(:,[1,1,2,2]);
tBM = sparse(brow,bcol, t_bM, size(elements.points,1),size(elements.points,1));

Z = sparse(row,col,Zsparse, size(elements.points,1),size(elements.points,1));
M = sparse(row,col,Msparse, size(elements.points,1),size(elements.points,1));
A = Z + 1i*beta*kappa*tBM;

b = tBM*hVec - M*FVec;


%% solve the system
U = A\b;

end

