function [modDomain, modBoundary, obs] = forwardOperatorWestervelt(elements, boundary, observation, timeMeshh, L, M, u, eta, b, s, gamma)

% the solution u is of the format time x space(evaluated at vertices of
% triangles)

% the easy one first, the observation
obs = u(:,observation);

%(b u - \eta u^2)_{tt}
uh = (b.'.*u - eta.'.*u.^2);
uh_t = (circshift(uh,-1,1) - uh)./timeMeshh;
uh_tt = (circshift(uh_t,-1,1)- uh_t)./timeMeshh;

Deltau = (M\L * u.').';
Deltau_t = (circshift(Deltau,-1,1) - Deltau)./timeMeshh;

% only in \Omega
fullboundaryIdx = elements.edges(:,1);
interiorIdx = setdiff(1:(size(elements.points,1)), fullboundaryIdx);

modfull = uh_tt - s.' .* Deltau - Deltau_t;
modDomain = modfull(:,interiorIdx);

end