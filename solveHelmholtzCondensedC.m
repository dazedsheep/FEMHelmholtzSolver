function [U] = solveHelmholtzCondensedC(elements, omega, gamma, kappa, beta, f, hI, n, K, rowK, colK, M_t, tBM)

hVec = zeros(n,1);

KappaSq = repmat(mean(kappa(elements.nodeIndex).^2,2),1,9,1).';
MC = M_t .* KappaSq;

Z_t = reshape(K, 9, size(K,3)) - MC;
% put the elements in Z_t on the right position with the local to global
% index
Z = sparse(rowK, colK, Z_t, size(elements.points,1),size(elements.points,1));

% sparse mass matrix
M = sparse(rowK,colK, M_t, size(elements.points,1),size(elements.points,1));

A = Z + (1i.*beta.*omega + gamma).*tBM;

% get the neumann/robin boundary values
%hVec(elements.bedges(:,1)) = hI(elements.bedges(:,1));
%f(elements.bedges(:,1)) = 0;
% right hand side
b = tBM*hVec + M*f;

%% solve the system
U = A\b;

end

