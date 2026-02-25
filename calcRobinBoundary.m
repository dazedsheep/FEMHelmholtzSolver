function [robinBoundary] = calcRobinBoundary( ...
    elements, u, gamma, Gx, Gy)

Nt = size(u,1);
Np = size(u,2);


% compute gradient for ALL time steps at once
Ux = (Gx * u.').';
Uy = (Gy * u.').';

nb = elements.boundaryIdx;

normal_x = elements.boundaryNormals(:,1).';
normal_y = elements.boundaryNormals(:,2).';

normalGrad = Ux(:,nb).*normal_x + Uy(:,nb).*normal_y;

robinBoundary = zeros(Nt,Np);
robinBoundary(:,nb) = gamma.*u(:,nb) + normalGrad;

end
