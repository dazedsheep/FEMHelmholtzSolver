function [modDomain, modBoundary, obs] = linearisedForwardOperatorAllAtOnce(elements, observation, timeMeshh, A, x0, dx, Gx, Gy)
Nt = size(x0.u0,1);
Np = size(x0.u0,2);
% this also support time dependent coefficients

uh = x0.b0 .* dx.du - 2.* x0.eta0 .* x0.u0 .* dx.du - dx.deta .* (x0.u0.^2) + dx.db.*x0.u0;
dutt = dtt(uh, timeMeshh);

s0laplace = x0.s0 .* laplacian(dx.du, A);
dut = dt(dx.du, timeMeshh);
dutlaplace = laplacian(dut, A);
dslaplace = dx.ds .* laplacian(x0.u0, A);

modDomain = dutt - s0laplace - dutlaplace - dslaplace;

%% Observation
obs = dx.du(:,observation);

% compute gradient for ALL time steps at once
Ux = (Gx * dx.du.').';
Uy = (Gy * dx.du.').';

nb = elements.boundaryIdx;

normal_x = elements.boundaryNormals(:,1).';
normal_y = elements.boundaryNormals(:,2).';

normalGrad = Ux(:,nb).*normal_x + Uy(:,nb).*normal_y;

modBoundary = zeros(Nt,Np);
modBoundary(:,nb) = dx.gamma.*dx.du(:,nb) + normalGrad;

end

