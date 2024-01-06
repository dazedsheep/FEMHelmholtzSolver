function [U] = solveHelmholtzC(elements, omega, kappa, gamma, beta, f, hI, g, n)
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

for i=1:size(elements.triangles,2)
    % prepare the prescribed h on the boundary
    [K, MC, M, F] = assembleFEMC(elements.triangles{i}, kappa, f, quadratureParameters);
    
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
 
    Z = K - MC;
    Zsparse(localNodeIdx) = Z;
    Msparse(localNodeIdx) = M;
    
end

for i=1:size(elements.bedges, 1)
    hVec(elements.bedges(i,1)) = hI(elements.points(elements.bedges(i,1), 1), elements.points(elements.bedges(i,1), 2));
end


% process the transparent (Neumann or Robin) boundary conditions - this is
% already vectorized
e_Vec = elements.points(elements.bedges(:,1),:) - elements.points(elements.bedges(:,2),:);
e_len = sqrt(sum(e_Vec.^2,2));
t_bM = e_len'/6 .* [1;2;2;1]; % boundary element mass matrix
brow = elements.bedges(:,[1,2,1,2]);
bcol = elements.bedges(:,[1,1,2,2]);
tBM = sparse(brow,bcol, t_bM, size(elements.points,1),size(elements.points,1));

Z = sparse(row,col,Zsparse, size(elements.points,1),size(elements.points,1));
M = sparse(row,col,Msparse, size(elements.points,1),size(elements.points,1));

A = Z + (1i*beta*omega+gamma)*tBM;

b = tBM*hVec - M*FVec;

%% solve the system
U = A\b;

end

