function [modDomain, modBoundary, obs] = forwardOperatorWestervelt( ...
    elements, observation, timeMeshh, L, M, u, eta, b, s, gamma, Gx, Gy)

Nt = size(u,1);
Np = size(u,2);

%% Observation
obs = u(:,observation);

%% Time derivatives (vectorized, no circshift)
u_t = [ (u(2,:) - u(1,:))/timeMeshh; 
        (u(3:end,:) - u(1:end-2,:))/(2*timeMeshh); 
        (u(end,:) - u(end-1,:))/timeMeshh ];

u_tt = [ (u_t(2,:) - u_t(1,:))/timeMeshh; 
        (u_t(3:end,:) - u_t(1:end-2,:))/(2*timeMeshh); 
        (u_t(end,:) - u_t(end-1,:))/timeMeshh ];

uh   = b.'.*u - eta.'.*u.^2;
uh_t = [ (uh(2,:) - uh(1,:))/timeMeshh; 
        (uh(3:end,:) - uh(1:end-2,:))/(2*timeMeshh); 
        (uh(end,:) - uh(end-1,:))/timeMeshh ];

uh_tt = [(uh_t(2,:) - uh_t(1,:))/timeMeshh; 
        (uh_t(3:end,:) - uh_t(1:end-2,:))/(2*timeMeshh); 
        (uh_t(end,:) - uh_t(end-1,:))/timeMeshh ];
%% Laplacian (factorize once if reused often)
A = M\L;
Laplaceu   = (A * u.').';
Laplaceu_t = (A * u_t.').';

%% Domain residual
modDomain = zeros(Nt,Np);
modfull = uh_tt - s.' .* Laplaceu - Laplaceu_t;
modDomain(:,elements.interiorIdx) = modfull(:,elements.interiorIdx);

%% --------- FAST BOUNDARY GRADIENT ---------
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