function [modDomain, modBoundary, obs] = forwardOperatorAllAtOnce( ...
    elements, observation, timeMeshh, A, u, eta, b, s, gamma, Gx, Gy)

Nt = size(u,1);
Np = size(u,2);

%% Observation
obs = u(:,observation);

%% Time derivatives 
u_t = dt(u, timeMeshh);
uh   = (b.*u.' - eta.*(u.^2).').';
uh_tt = dtt(uh, timeMeshh);

%% Laplacians
sLaplaceu = s.'.* laplacian(u,A);
Laplaceu_t = laplacian(u_t, A);

%% Domain residual
modDomain = zeros(Nt,Np);
modfull = uh_tt - sLaplaceu - Laplaceu_t;
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