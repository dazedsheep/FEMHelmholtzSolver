function [U] = solveHelmholtzLM(elements, omega, gamma, kappa, beta, f, hI, n, K, rowK, colK, M_t, tBM)

hVec = hI;
KappaSq = repmat(mean(kappa(elements.nodeIndex),2),1,9,1).';
KappaSqSparse = sparse(rowK,colK, KappaSq, size(elements.points,1),size(elements.points,1));

MC = M_t .* KappaSqSparse;

A = K - MC + (1i.*beta.*omega + gamma).*tBM;

% get the neumann/robin boundary values
%hVec(elements.bedges(:,1)) = hI(elements.bedges(:,1));
%f(elements.bedges(:,1)) = 0;
% right hand side
b = tBM*hVec + M_t *f;

%% solve the system
U = A\b;

end

