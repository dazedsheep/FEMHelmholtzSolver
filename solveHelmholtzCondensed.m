function [U] = solveHelmholtzCondensed(elements, omega, gamma, kappa, beta, f, hI, n, K, rowK, colK, M_t, tBM)

hVec = zeros(n,1);
FVec = zeros(n,1);

Z_t = reshape(K, 9, size(K,3))  - kappa^2 * M_t;
% put the elements in Z_t on the right position with the local to global
% index
Z = sparse(rowK, colK, Z_t, size(elements.points,1),size(elements.points,1));

% sparse mass matrix
M = sparse(rowK,colK, M_t, size(elements.points,1),size(elements.points,1));

A = Z + (1i*beta*omega+gamma)*tBM;

% get the right hand side
FVec(unique(elements.nodeIndex(:,1)),1) = f(unique(elements.nodeIndex(:,1)));

% get the neumann/robin boundary values
hVec(elements.bedges(:,1)) = hI(elements.bedges(:,1));

% right hand side
b = tBM*hVec - M*FVec;

%% solve the system
U = gather(A\b);

end

