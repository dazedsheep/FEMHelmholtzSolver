function [modDomain, modBoundary, obs] = forwardOperatorAllAtOnce( ...
    elements, observation, timeMeshh, L, M, u, eta, b, s, gamma, Gx, Gy)

Nt = size(u,1);
Np = size(u,2);

%% Observation
obs = u(:,observation);

%% Time derivatives (vectorized)
u_t = zeros(Nt,Np);
u_t = [ (u(end,:) - u(1,:))/timeMeshh; 
        (u(3:end,:) - u(1:end-2,:))/(2*timeMeshh); 
        (u(1,:) - u(end,:))/timeMeshh ];
uh   = (b.*u.' - eta.*(u.^2).').';
uh_tt = [ (uh(3,:) - 2*uh(2,:) + uh(1,:)) / timeMeshh^2;
          (uh(3:end,:) - 2*uh(2:end-1,:) + uh(1:end-2,:)) / timeMeshh^2;
          (uh(end,:) - 2*uh(end-1,:) + uh(end-2,:)) / timeMeshh^2 ];

%% Laplacian using our stiffness matrix
A = M\L;
for i = 1: Nt
    sLaplaceu(i,:)   = s.*(A * u(i,:).');
    Laplaceu_t(i,:) = (A * u_t(i,:).');
end

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