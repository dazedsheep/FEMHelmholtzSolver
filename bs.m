%% residual test of forward operator
u0 = zeros(size(timeMesh,2), size(elements.points,1));
F = @(t,x,y) ...
    - b .* omega1.^2 .* (x.^2 + y.^2 + 1) .* cos(omega1 .* t) ...
    - 2 .* eta .* omega1.^2 .* (x.^2 + y.^2 + 1).^2 .* ...
      (1 - 2 .* cos(omega1 .* t).^2 - 2 .* cos(omega1 .* t)) ...
    - 4 .* s .* (cos(omega1 .* t) + 2) ...
    + 4 .* omega1 .* sin(omega1 .* t);
Fu  = @(t,x,y) ...
   omega1.^2.*(x.^2 + y.^2 + 1).*(cos(omega1.*t).*(4.*eta.*(x.^2 + y.^2 + 1) - b) + 2.*eta.*(x.^2 + y.^2 + 1).* cos(2.*t.*omega1));

Flapu = @(t,x,y) (4 .* (cos(omega1 .* t) + 2));
Flaput = @(t,x,y) (- 4 .* omega1 .* sin(omega1 .* t));
Fu_t = @(t,x,y)(-(x.^2 + y.^2 + 1).*sin(omega1.*t).*omega1);


analyticModDomain = zeros(size(timeMesh,2), size(elements.points,1));

for i=1:size(timeMesh,2)
    u0(i,:) = u1(timeMesh(i),elements.points(:,1), elements.points(:,2));
    uh(i,:) = (b.*u1(timeMesh(i),elements.points(:,1), elements.points(:,2)) - eta.*(u1(timeMesh(i),elements.points(:,1), elements.points(:,2)).^2) );
    auh(i,:) = Fu(timeMesh(i),elements.points(:,1), elements.points(:,2));
    alapu(i,:) = Flapu(timeMesh(i),elements.points(:,1), elements.points(:,2));
    alaput(i,:) = Flaput(timeMesh(i),elements.points(:,1), elements.points(:,2));
    u_t(i,:) = Fu_t(timeMesh(i),elements.points(:,1), elements.points(:,2));
    analyticModDomain(i,:) = F(timeMesh(i),elements.points(:,1), elements.points(:,2));
end

%%
% ---------------------------------------------------------
% Parameters
% ---------------------------------------------------------
R      = 0.2;

Nx = 201;          % spatial resolution
Nt = length(timeMesh);
% ---------------------------------------------------------
% Spatial grid
% ---------------------------------------------------------
x = linspace(-R,R,Nx);
y = linspace(-R,R,Nx);
[X,Y] = meshgrid(x,y);

h = x(2) - x(1);

mask = (X.^2 + Y.^2) <= R^2;

% ---------------------------------------------------------
% Allocate space-time array
% u is Nx x Nx x Nt
% ---------------------------------------------------------
u = zeros(Nx,Nx,Nt);
A = X.^2 + Y.^2 + 1;

% ---------------------------------------------------------
% Evaluate u on space-time mesh
% ---------------------------------------------------------
for k = 1:Nt
    u(:,:,k) = A .* (cos(omega1*timeMesh(k)) + 2);
end

% ---------------------------------------------------------
% Compute Laplacian for each time
% ---------------------------------------------------------
Lap_u = zeros(size(u));

for k = 1:Nt
    
    uk = u(:,:,k);
    
    uxx = zeros(Nx,Nx);
    uyy = zeros(Nx,Nx);
    
    % central differences
    uxx(:,2:end-1) = (uk(:,3:end) - 2*uk(:,2:end-1) + uk(:,1:end-2)) / h^2;
    uyy(2:end-1,:) = (uk(3:end,:) - 2*uk(2:end-1,:) + uk(1:end-2,:)) / h^2;
    
    Lap_u(:,:,k) = uxx + uyy;
end

% Apply disk mask
for k = 1:Nt
    tmp = Lap_u(:,:,k);
    tmp(~mask) = NaN;
    Lap_u(:,:,k) = tmp;
end

% ---------------------------------------------------------
% Exact Laplacian (for verification)
% ---------------------------------------------------------
Lap_exact = zeros(size(u));

for k = 1:Nt
    Lap_exact(:,:,k) = 4*(cos(omega1*timeMesh(k))+2);
end

Lap_exact(~repmat(mask,1,1,Nt)) = NaN;

% ---------------------------------------------------------
% Error
% ---------------------------------------------------------
err = max(abs(Lap_u(:) - Lap_exact(:)),[],'omitnan');
fprintf('Max space-time Laplacian error: %.3e\n',err);


%%
%  
A = M\L;
[modDomain, modBoundary, obs] = forwardOperatorAllAtOnce(elements, measurementPointsIdx, timeMeshh, A, u0, eta, b, s, gamma, Gx, Gy);

[modDomainSol, modBoundarySol, obsSol] = forwardOperatorAllAtOnce(elements, measurementPointsIdx, timeMeshh, A, u1s, eta, b, s, gamma, Gx, Gy);

l2boundaryValError = sqrt(sum(sum(abs(modBoundary - modBoundarySol).^2,2),1));
linfboundaryValError = sqrt(max(max(abs(modBoundary - modBoundarySol))));
%% TODO: implement linearised forward operator and its adjoint + min prob

% setup initial state(s)

% sample the reference states on our triangle mesh

for i=1:size(timeMesh,2)
    u0_1(i,:) = u1(timeMesh(i), elements.points(:,1), elements.points(:,2));
    u0_2(i,:) = u2(timeMesh(i), elements.points(:,1), elements.points(:,2));
    u0_3(i,:) = u3(timeMesh(i), elements.points(:,1), elements.points(:,2));
end

% assemble the measurement vector
[modDomainSol1, modBoundarySol1, obsSol1] = forwardOperatorAllAtOnce(elements, measurementPointsIdx, timeMeshh, A, u1s, eta, b, s, gamma, Gx, Gy);
[modDomainSol2, modBoundarySol2, obsSol2] = forwardOperatorAllAtOnce(elements, measurementPointsIdx, timeMeshh, A, u2s, eta, b, s, gamma, Gx, Gy);
[modDomainSol3, modBoundarySol3, obsSol3] = forwardOperatorAllAtOnce(elements, measurementPointsIdx, timeMeshh, A, u3s, eta, b, s, gamma, Gx, Gy);

%%
% in what point do we linearise (we will also start here)
s0 = ones(size(elements.points(:,1))).*500;
b0 = ones(size(elements.points(:,1)));
eta0 = 0;

% for each outer iteration we have to solve a system Ax = b, which we can't
% solve directly --> min prob in inner loop
x0.s0 = s0.';
x0.b0 = b0.';
x0.eta0 = 0;
x0.gamma = gamma;

x_n.u1 = u0_1;
x_n.u2 = u0_2;
x_n.u3 = u0_3;

x_n.u1 = u1s;
x_n.u2 = u2s;
x_n.u3 = u3s;

x_n.eta = x0.eta0.';
x_n.b = x0.b0.';
x_n.s = x0.s0.';
x_n.gamma = gamma;

alpha_n = 1;

% we need to compute b_n = K*(h - F(x_n)) - P*P(x_n) + \alpha_n (x_0 -
% x_n)
[modDomain1, ~, obs1] = forwardOperatorAllAtOnce(elements, measurementPointsIdx, timeMeshh, A, x_n.u1, x_n.eta, x_n.b, x_n.s, x_n.gamma, Gx, Gy);
[modDomain2, ~, obs2] = forwardOperatorAllAtOnce(elements, measurementPointsIdx, timeMeshh, A, x_n.u2, x_n.eta, x_n.b, x_n.s, x_n.gamma, Gx, Gy);
[modDomain3, ~, obs3] = forwardOperatorAllAtOnce(elements, measurementPointsIdx, timeMeshh, A, x_n.u3, x_n.eta, x_n.b, x_n.s, x_n.gamma, Gx, Gy);

hx.modDomain = modDomainSol1 - modDomain1;
hx.obs = measurement_u1 - obs1;
x0.u0 = u0_1;
[adju1] = adjointLinearisedForwardOperatorAllAtOnce(elements, timeMeshh, A, x0, hx, Gx, Gy);
hx.modDomain = modDomainSol2 - modDomain2;
hx.obs = measurement_u2 - obs2;
x0.u0 = u0_2;
[adju2] = adjointLinearisedForwardOperatorAllAtOnce(elements, timeMeshh, A, x0, hx, Gx, Gy);
hx.modDomain = modDomainSol3 - modDomain3;
hx.obs = measurement_u3 - obs3;
x0.u0 = u0_3;
[adju3] = adjointLinearisedForwardOperatorAllAtOnce(elements, timeMeshh, A, x0, hx, Gx, Gy);

adju1 = projectTimeConstant(timeMeshh, adju1);
adju2 = projectTimeConstant(timeMeshh, adju2);
adju3 = projectTimeConstant(timeMeshh, adju3);

b_n.u1 = adju1.du - x_n.u1 + alpha_n .* (u0_1 - x_n.u1);
b_n.u2 = adju2.du - x_n.u2 + alpha_n .* (u0_2 - x_n.u2);
b_n.u3 = adju3.du - x_n.u3 + alpha_n .* (u0_3 - x_n.u3);
b_n.eta = 1/3 .* (adju1.deta + adju2.deta + adju2.deta);
b_n.b = 1/3 .* (adju1.db + adju2.db + adju2.db);
b_n.s = 1/3 .* (adju1.ds + adju2.ds + adju2.ds);
b_n.obs1 = adju1.dobs - x_n.u1(:,measurementPointsIdx) + alpha_n .*(u0_1(:,measurementPointsIdx) - x_n.u1(:,measurementPointsIdx));
b_n.obs2 = adju2.dobs - x_n.u2(:,measurementPointsIdx) + alpha_n .*(u0_2(:,measurementPointsIdx) - x_n.u2(:,measurementPointsIdx));
b_n.obs3 = adju3.dobs - x_n.u3(:,measurementPointsIdx) + alpha_n .*(u0_3(:,measurementPointsIdx) - x_n.u3(:,measurementPointsIdx));
tau = 10e-6;

% (K*K + P*P + alpha_n)(d_n)

d_n.u1 = zeros(size(u1s));
d_n.u2 = zeros(size(u1s));
d_n.u3 = zeros(size(u1s));
d_n.obs1 = zeros(size(obsSol1));
d_n.obs2 = zeros(size(obsSol1));
d_n.obs3 = zeros(size(obsSol1));
d_n.eta = 0;
d_n.s = 0;
d_n.b = 0;
for i=1:20
% A = K*K + P*P + alpha_n
dx.du = d_n.u1;
dx.ds = d_n.s;
dx.db = d_n.b;
dx.deta = d_n.eta;
dx.gamma = gamma;
x0.u0 = u0_1;
[modDomain, modBoundary, obs] = linearisedForwardOperatorAllAtOnce(elements, measurementPointsIdx, timeMeshh, A, x0, dx, Gx, Gy);
hx.modDomain = modDomain;
hx.obs = obs;
[adju1] = adjointLinearisedForwardOperatorAllAtOnce(elements, timeMeshh, A, x0, hx, Gx, Gy);
dx.du = d_n.u2;
x0.u0 = u0_2;
[modDomain, modBoundary, obs] = linearisedForwardOperatorAllAtOnce(elements, measurementPointsIdx, timeMeshh, A, x0, dx, Gx, Gy);
hx.modDomain = modDomain;
hx.obs = obs;
[adju2] = adjointLinearisedForwardOperatorAllAtOnce(elements, timeMeshh, A, x0, hx, Gx, Gy);
dx.du = d_n.u3;
x0.u0 = u0_3;
[modDomain, modBoundary, obs] = linearisedForwardOperatorAllAtOnce(elements, measurementPointsIdx, timeMeshh, A, x0, dx, Gx, Gy);
hx.modDomain = modDomain;
hx.obs = obs;
[adju3] = adjointLinearisedForwardOperatorAllAtOnce(elements, timeMeshh, A, x0, hx, Gx, Gy);

adju1 = projectTimeConstant(timeMeshh, adju1);
adju2 = projectTimeConstant(timeMeshh, adju2);
adju3 = projectTimeConstant(timeMeshh, adju3);
% b_n - Ad_n
deta = b_n.eta - (1/3 .* (adju1.deta + adju2.deta + adju2.deta) + alpha_n .* d_n.eta);
db = b_n.b - (1/3 .* (adju1.db + adju2.db + adju2.db) + alpha_n .* d_n.b);
ds = b_n.s - (1/3 .* (adju1.ds + adju2.ds + adju2.ds) + alpha_n .* d_n.s);
du1 = b_n.u1 - (adju1.du + alpha_n .* d_n.u1);
du2 = b_n.u2 - (adju2.du + alpha_n .* d_n.u2);
du3 = b_n.u3 - (adju3.du + alpha_n .* d_n.u3);
dobs1 = b_n.obs1 - (adju1.dobs + alpha_n .* d_n.obs1);
dobs2 = b_n.obs2 - (adju1.dobs + alpha_n .* d_n.obs1);
dobs3 = b_n.obs3 - (adju1.dobs + alpha_n .* d_n.obs1);

d_n.u1 = d_n.u1 - tau.*du1;
d_n.u2 = d_n.u2 - tau.*du2;
d_n.u3 = d_n.u3 - tau.*du3;
d_n.obs1 = d_n.obs1 - tau.*dobs1;
d_n.obs2 = d_n.obs2 - tau.*dobs2;
d_n.obs3 = d_n.obs3 - tau.*dobs3;
d_n.eta = d_n.eta - tau.*deta;
d_n.s = d_n.s - tau.*ds;
d_n.b = d_n.b - tau.*db;
end