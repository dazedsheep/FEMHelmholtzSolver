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
timeMesh = linspace(0,1/f1,1/timeMeshh);

% our domain
bcenter = [0,0];
brad = 0.2;
domain = [bcenter, brad];

% specify the mesh parameter
meshSize = 0.005;

% compute the triangle mesh
[elements] = initializeMultiLeveLSolver(meshSize, domain);

% pre-compute some useful things w.r.t. the triangular mesh
fullboundaryIdx = elements.edges(:,1);
interiorIdx = setdiff(1:(size(elements.points,1)), fullboundaryIdx);
elements.interiorIdx = interiorIdx;
elements.boundaryIdx = fullboundaryIdx;
elements.boundaryNormals = 1./sqrt(sum(elements.points(fullboundaryIdx,:).^2,2)).*elements.points(fullboundaryIdx,:); % our center is (0,0), so -> normalisation is suffices


% specify the parameters we want to reconstruct
% boundary parameters (not reconstructed)
gamma = 1;

beta = 0;   % this is important check paper for clarification

% define a phantom in our domain with different speed of sound, diffusivity
% and nonlinearity parameter
diffusivity = 0.5;
values = [5]; % B/A of phantoms
radii = [0.05];
diffusivityPhantoms = [5]; % this allows to adjust the diffusivity for the phantoms
centers = [0; 0];

massDensity = 1000; %kg/m^3

speed_of_sound = 1;

N = 10; % number of harmonics-1 we will compute

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
[cN, U1, F1] = solveWesterveltMultiLevelBoundaryExcitation(elements, omega1, beta, gamma, squeeze(kappasq(:,:,1)), squeeze(excitations(:,:,1)), eta, b, 15, N, 10^(-12));
[cN, U2, F2] = solveWesterveltMultiLevelBoundaryExcitation(elements, omega2, beta, gamma, squeeze(kappasq(:,:,2)), squeeze(excitations(:,:,2)), eta, b, 15, N, 10^(-12));
[cN, U3, F3] = solveWesterveltMultiLevelBoundaryExcitation(elements, omega1, beta, gamma, squeeze(kappasq(:,:,3)), squeeze(excitations(:,:,3)), eta, b, 15, N, 10^(-12));
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
L = cotmatrix(elements.points, elements.tri); % laplacian matrix (space)
M = massmatrix(elements.points, elements.tri); % mass matrix (space)

[modDomain, modBoundary, obs] = forwardOperatorWestervelt(elements, measurementPointsIdx, timeMeshh, L, M, u1s, eta, b, s, gamma);

%testu1boundary = zeros(1,size(u2s,2));
% testu1boundary(1,boundaryPointsSourceIdx) = gamma.*u1(0, boundaryPointsSource(:,1), boundaryPointsSource(:,2)) + dot(squeeze(u1grad(0,boundaryPointsSource(:,1),boundaryPointsSource(:,2))).', boundaryPointsSourceNormals.').';