function [modDomain, modBoundary, obs] = forwardOperatorAllAtOnce( ...
    elements, observation, timeMeshh, L, M, u, eta, b, s, gamma, Gx, Gy)

Nt = size(u,1);
Np = size(u,2);
% local to global indices
n = size(elements.points,1);
rowK = elements.tri(:, [1 2 3 1 2 3 1 2 3]).';
colK = elements.tri(:, [1 1 1 2 2 2 3 3 3]).';
etaFEM = repmat(mean(eta(elements.nodeIndex),2),1,9,1).';
bFEM = repmat(mean(b(elements.nodeIndex),2),1,9,1).';
sFEM = repmat(mean(s(elements.nodeIndex),2),1,9,1).';
K = sparse(rowK, colK, L, n, n);
Meta = M .* etaFEM;
Mb = M .* bFEM;
Ks = L .* sFEM;
Meta = sparse(rowK, colK, Meta, n,n);
Mb = sparse(rowK, colK, Mb, n,n);
Ks = sparse(rowK, colK, Ks, n,n);
%% Observation
obs = u(:,observation);

%% Time derivatives (vectorized, no circshift)
u_t = [ (u(2,:) - u(1,:))/timeMeshh; 
        (u(3:end,:) - u(1:end-2,:))/(2*timeMeshh); 
        (u(end,:) - u(end-1,:))/timeMeshh ];
uh   = (Mb*u.' - Meta*(u.^2).').';
uh_t = [ (uh(2,:) - uh(1,:))/timeMeshh; 
        (uh(3:end,:) - uh(1:end-2,:))/(2*timeMeshh); 
        (uh(end,:) - uh(end-1,:))/timeMeshh ];
uh_tt = [(uh_t(2,:) - uh_t(1,:))/timeMeshh; 
        (uh_t(3:end,:) - uh_t(1:end-2,:))/(2*timeMeshh); 
        (uh_t(end,:) - uh_t(end-1,:))/timeMeshh ];
%% Laplacian (factorize once if reused often)
%A = M\L;
Laplaceu   = (Ks * u.').';
Laplaceu_t = (K * u_t.').';
%% Domain residual
modDomain = zeros(Nt,Np);
modfull = uh_tt -  Laplaceu - Laplaceu_t;
modDomain(:,elements.interiorIdx) = modfull(:,elements.interiorIdx);

% compute gradient for ALL time steps at once

Ux = (Gx * u.').';
Uy = (Gy * u.').';

nb = elements.boundaryIdx;

normal_x = elements.boundaryNormals(:,1).';
normal_y = elements.boundaryNormals(:,2).';

normalGrad = Ux(:,nb).*normal_x + Uy(:,nb).*normal_y;

modBoundary = zeros(Nt,Np);
modBoundary(:,nb) = gamma.*u(:,nb) + normalGrad;

end