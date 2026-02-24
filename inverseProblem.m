clear all

% specify our reference states
f1 = 10;    % 10 Hz
f2 = 20;    % 20 Hz             
omega1 = 2*pi*f1;     
omega2 = 2*pi*f2;     

u1 = @(t,x,y) (x.^2 + y.^2 + 1) .* (cos(omega1 .* t) + 2);
u2 = @(t,x,y) (x.^2 + y.^2 + 1) .* (cos(omega2 .* t) + 2);
u3 = @(t,x,y) 2 .* u1(t,x,y);
u1grad = @(t,x,y)  cat(3,...
    (2.*x).* (cos(omega1 .* t) + 2), ...
    (2.*y).* (cos(omega1 .* t) + 2));
u1F = @(x,y) (x.^2 + y.^2 + 1); 
u1C = @(x,y) (x.^2 + y.^2 + 1) .* 2;
u1Fgrad = @(x,y) cat(3,...
    2 .* x, ...
    2 .* y);
u1Cgrad = @(x,y) 2.*u1Fgrad(x,y);

% specify our time space cylinder and calculate the triangle mesh in space
% and the mesh in time
timeMeshh = 0.01;
% time domain (lowest frequency determines the duration)
timeMesh = linspace(0,1/f1,1/timeMeshh - 1);
% recompute time difference
timeMeshh = timeMesh(2) - timeMesh(1); % careful this is the time diff!

% our domain
bcenter = [0,0];
brad = 0.2;
domain = [bcenter, brad];

% specify the mesh parameter
meshSize = 0.01;

% compute the triangle mesh
[elements] = initializeMultiLeveLSolver(meshSize, domain);

% pre-compute some useful things w.r.t. the triangular mesh
fullboundaryIdx = elements.edges(:,1);
interiorIdx = setdiff(1:(size(elements.points,1)), fullboundaryIdx);
elements.interiorIdx = interiorIdx;
elements.boundaryIdx = fullboundaryIdx;
elements.boundaryNormals = 1./sqrt(sum(elements.points(fullboundaryIdx,:).^2,2)).*elements.points(fullboundaryIdx,:); % our center is (0,0), so -> normalisation suffices


% specify the parameters we want to reconstruct
% boundary parameters (not reconstructed)
gamma = 1;

beta = 0;   % this is important check paper for clarification

% define a phantom in our domain with different speed of sound, diffusivity
% and nonlinearity parameter
diffusivity = 0.005;
values = [5]; % B/A of phantoms
radii = [0.05];
diffusivityPhantoms = [20]; % this allows to adjust the diffusivity for the phantoms
centers = [0; 0];

massDensity = 1000; %kg/m^3

speed_of_sound = 2;

N = 5; % number of harmonics-1 we will compute

% create the space dependent parameters
sourceValueDomain = 2; % B/A of domain

eta = constructNonlinearityDivB(elements, massDensity, speed_of_sound, diffusivity, diffusivityPhantoms, centers, radii, values, sourceValueDomain, true); %nonlinearity scaled by 1/b
s = constructSquaredSpeedOfSoundDivB(elements, speed_of_sound, diffusivity, diffusivityPhantoms, centers, radii); % speed of sound scaled by 1/b

b = constructReciprocalDiffusivity(elements, diffusivity, diffusivityPhantoms, centers, radii);

% the complex wavenumber, here we compute the square wave number 
kappasq = constructKappaReparameterized(elements, s, b, [omega1 omega2 omega1], N); % compute all the complex wave numbers needed

% since u^0_j  does not need to solve the PDE,... 
% construct the boundary excitations (same source just with different
% frequency

%% prepare the source(s)
sourceEdge = 1; % we impose the source on the boundary (negative quadrant)
%boundaryPointsSourceIdx = elements.edges((elements.edges(:,3) == sourceEdge),1);

boundaryPointsSourceIdx = elements.boundaryIdx;

boundaryPointsSource = elements.points(boundaryPointsSourceIdx,:);

% compute the normals
boundaryPointsSourceNormals = 1./sqrt(sum(elements.points(boundaryPointsSourceIdx,:).^2,2)).*elements.points(boundaryPointsSourceIdx,:); % our center is (0,0), so -> normalisation is suffices

% frequency part (robin boundary)
sourceFrequency = zeros(size(elements.points,1),1);
sourceFrequency(boundaryPointsSourceIdx) = 1/2.*(gamma.*u1F(boundaryPointsSource(:,1),boundaryPointsSource(:,2)) + dot(squeeze(u1Fgrad(boundaryPointsSource(:,1),boundaryPointsSource(:,2))).',boundaryPointsSourceNormals.').');

% constant part (robin boundary)
sourceConstant = zeros(size(elements.points,1),1);
sourceConstant(boundaryPointsSourceIdx) = gamma.*u1C(boundaryPointsSource(:,1), boundaryPointsSource(:,2)) + dot(squeeze(u1Cgrad(boundaryPointsSource(:,1),boundaryPointsSource(:,2))).',boundaryPointsSourceNormals.').';

%%
excitations = zeros(size(elements.points,1), N, 3);
excitations(:,1,1) = sourceConstant;
excitations(:,2,1) = sourceFrequency;
excitations(:,1,2) = sourceConstant;
excitations(:,2,2) = sourceFrequency;
excitations(:,1,3) = 2.*sourceConstant;
excitations(:,2,3) = 2.*sourceFrequency;
%%
[cN, U1, F1] = solveWesterveltMultiLevelBoundaryExcitation(elements, omega1, beta, gamma, squeeze(kappasq(:,:,1)), squeeze(excitations(:,:,1)), eta, b, 5, N, 10^(-12));
[cN, U2, F2] = solveWesterveltMultiLevelBoundaryExcitation(elements, omega2, beta, gamma, squeeze(kappasq(:,:,2)), squeeze(excitations(:,:,2)), eta, b, 5, N, 10^(-12));
[cN, U3, F3] = solveWesterveltMultiLevelBoundaryExcitation(elements, omega1, beta, gamma, squeeze(kappasq(:,:,3)), squeeze(excitations(:,:,3)), eta, b, 5, N, 10^(-12));
%%
% compute the solutions on the time - space mesh

u1s = calcSolution(timeMesh, squeeze(U1(N,:,:)), omega1);
u2s = calcSolution(timeMesh, squeeze(U2(N,:,:)), omega2);
u3s = calcSolution(timeMesh, squeeze(U3(N,:,:)), omega1);

%%
% define \Sigma our measurement manifold/discrete points
% the triangulation already defines the the edges of our doimain 1:4
% (circle sectors)

measurementEdge = 3; % positive quadrant edge

% fetch the boundary points
measurementPointsIdx = elements.edges((elements.edges(:,3) == measurementEdge),1);

measurement_u1 = u1s(:,measurementPointsIdx);
measurement_u2 = u2s(:,measurementPointsIdx);
measurement_u3 = u3s(:,measurementPointsIdx);
%%
% There is a bunch of things we can prepare beforehand computation
L = cotmatrix(elements.points, elements.tri); % laplacian matrix (space) [stiffnes matrix]
M = massmatrix(elements.points, elements.tri); % mass matrix (space)

% pre compoute the gradient operator
p = elements.points;
t = elements.tri(:,1:3);

x = p(:,1);
y = p(:,2);

nt = size(t,1);
np = size(p,1);

% triangle vertex indices
i1 = t(:,1); i2 = t(:,2); i3 = t(:,3);

% coordinates
x1 = x(i1); y1 = y(i1);
x2 = x(i2); y2 = y(i2);
x3 = x(i3); y3 = y(i3);

% triangle areas
area = 0.5*((x2-x1).*(y3-y1)-(x3-x1).*(y2-y1));

% gradients of basis functions
b1 = (y2-y3)./(2*area);
b2 = (y3-y1)./(2*area);
b3 = (y1-y2)./(2*area);

c1 = (x3-x2)./(2*area);
c2 = (x1-x3)./(2*area);
c3 = (x2-x1)./(2*area);

% assemble sparse gradient matrices
rows = repmat((1:nt)',1,3);

Gx_elem = sparse(rows(:), t(:), [b1;b2;b3], nt, np);
Gy_elem = sparse(rows(:), t(:), [c1;c2;c3], nt, np);

% average to nodes (mass lumping style)
Mmap = sparse(repmat((1:nt)',3,1), t(:), 1, nt, np);
weights = Mmap' * abs(area);

Gx = (Mmap' * (abs(area).*Gx_elem)) ./ weights;
Gy = (Mmap' * (abs(area).*Gy_elem)) ./ weights;

%testu1boundary = zeros(1,size(u2s,2));
% testu1boundary(1,boundaryPointsSourceIdx) = gamma.*u1(0, boundaryPointsSource(:,1), boundaryPointsSource(:,2)) + dot(squeeze(u1grad(0,boundaryPointsSource(:,1),boundaryPointsSource(:,2))).', boundaryPointsSourceNormals.').';
% this is the check for the scaled version:
%figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2), real(-squeeze(kappasq(:,2,1)).*squeeze(U1(N,2,:)) - M\L * squeeze(U1(N,2,:))), 'facecolor', 'interp'); shading interp;
% this is the check for the unscaled version:
%figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2), real(-squeeze(-b.*squeeze(U1(N,2,:))) - (s + 1i .*1.*omega1) .* M\L * squeeze(U1(N,2,:))), 'facecolor', 'interp'); shading interp;

%% sanity check
% bn = elements.boundaryIdx;
% harmonicIdx = 2;
% Ux = (Gx * squeeze(U1(N,harmonicIdx,:))).';  % N x 1
% Uy = (Gy * squeeze(U1(N,harmonicIdx,:))).';
% normal_x = elements.boundaryNormals(:,1);
% normal_y = elements.boundaryNormals(:,2);
% gradNormal =  (Ux(bn).' .* normal_x + Uy(bn).' .* normal_y);
% excitationCheck = zeros(size(elements.points,1),1);
% excitationCheck(bn) = gamma*squeeze(U1(N,harmonicIdx,bn)) + gradNormal;
% % check the domain residual, expect the biggest error on the boundary
% domresidual = squeeze(-Mass*kappasq(:,harmonicIdx,1)).*squeeze(U1(N,harmonicIdx,:)) + (S * squeeze(U1(N,harmonicIdx,:))) + MassB * squeeze(U1(N,harmonicIdx,:)) - Mass * F1(harmonicIdx,:).';
% u1boundary(elements.boundaryIdx) = gamma.*u1(0,boundaryPointsSource(:,1), boundaryPointsSource(:,2)) + dot(squeeze(u1grad(0,boundaryPointsSource(:,1),boundaryPointsSource(:,2))).',boundaryPointsSourceNormals.').';
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
%% TODO: implement linearised forward operator and its adjoint
x0.u0 = u0;
x0.s0 = s.';
x0.b0 = b.';
x0.eta0 = 0;
x0.gamma = gamma;

dx.du = u1s - u0;
dx.ds = 0;
dx.db = 0;
dx.deta = 0;
dx.gamma = gamma;

[modDomainK, modBoundaryK, obsK] = linearisedForwardOperatorAllAtOnce(elements, measurementPointsIdx, timeMeshh, A, x0, dx, Gx, Gy);

%%
point = [0.0;0.05];

[v,idx] = min(sum((elements.points - point(:)').^2,2)); 

node = [elements.points(idx,1);elements.points(idx,2)];

T = 1/f2;
% sampling frequency in time
Fs = 1/T * 2 * (N);
omega = omega1;
U = squeeze(U1(N,:,:));

Ns = 2000;
pC = zeros(1,Ns);
for m=0:(N-1)
    pC = pC + U(m+1,idx) .* exp(1i.*m.*omega.*(0:(Ns-1))*1/Fs);
end


MaxBins = 5;
P0 = max(max(u1s(:,:)));

window = hanning(Ns);
freq =  Fs/Ns*(0:(Ns/2));
y = abs(fft(window'.*real(pC)))/Ns;
y1 = y(1:Ns/2+1);
y1(2:end-1) = 2*y1(2:end-1);
shiftedTF = fftshift(fft(real(pC)))/Ns;
TFdB = 10*log10(y1/P0);
fscaling = 10^3;
M = min(MaxBins, Ns/2 + 1);
xaxis = freq./fscaling;
figure, plot(xaxis, TFdB(1:Ns/2+1))
xlabel("Frequency [kHz]")
ylabel("P/P0 [dB]")
