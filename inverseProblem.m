clear all

% specify our reference states
f1 = 100;    % 10 Hz
f2 = 200;    % 20 Hz             
omega1 = 2*pi*f1;     
omega2 = 2*pi*f2;     

u1 = @(t,x,y) (x.^2 + y.^2 + 1) .* (cos(omega1 .* t) + 2);
u2 = @(t,x,y) (x.^2 + y.^2 + 1) .* (cos(omega2 .* t) + 2);
u3 = @(t,x,y) 2 .* u1(t,x,y);

u1F = @(x,y) (x.^2 + y.^2 + 1); 
u1C = @(x,y) (x.^2 + y.^2 + 1) .* 2;
u1Fgrad = @(x,y) cat(3,...
    2 .* x, ...
    2 .* y);
u1Cgrad = @(x,y) 2.*u1Fgrad(x,y);

u1grad = @(t,x,y) cat(3, ...
    2 .* x .* (cos(omega1 .* t) + 2), ...
    2 .* y .* (cos(omega1 .* t) + 2));
u2grad = @(t,x,y) cat(3, ...
    2 .* x .* (cos(omega2 .* t) + 2), ...
    2 .* y .* (cos(omega2 .* t) + 2));
u3grad = @(t,x,y) 2.* u1grad(t,x,y);

% specify our time space cylinder and calculate the triangle mesh in space
% and the mesh in time

% time domain (lowest frequency determines the duration)
timeMesh = linspace(0,1/f1,100);

% our domain
bcenter = [0,0];
brad = 0.2;
domain = [bcenter, brad];

% specify the mesh parameter
meshSize = 0.01;

% compute the triangle mesh
[elements] = initializeMultiLeveLSolver(meshSize, domain);

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

speed_of_sound = 1480;

N = 10; % number of harmonics we will compute

% create the space dependent parameters
sourceValueDomain = 2; % B/A of domain

eta = constructNonlinearityDivB(elements, massDensity, speed_of_sound, diffusivity, diffusivityPhantoms, centers, radii, values, sourceValueDomain, true); %nonlinearity scaled by 1/b

s = constructSquaredSpeedOfSoundDivB(elements, speed_of_sound, diffusivity, diffusivityPhantoms, centers, radii); % speed of sound scaled by 1/b

b = constructReciprocalDiffusivity(elements, diffusivity, diffusivityPhantoms, centers, radii);

% the complex wavenumber 
kappa = constructKappaReparameterized(elements, s, b, [omega1 omega2 omega1], N); % compute all the complex wave numbers needed

% since u^0_j  does not need to solve the PDE,... 
% construct the boundary excitations (same source just with different
% frequency

% realistic piezoelectric elements are not of infinitesimal small size
excitationPoints = [0.0,0.0];
pressure = 10000;
excitationPointsSize = [0.001];
excitations = zeros(size(elements.points,1), N, 3);

%% prepare the source(s)
sourceEdge = 1; % we impose the source on the boundary (negative quadrant)
boundaryPointsSourceIdx = elements.edges((elements.edges(:,3) == sourceEdge),1);

boundaryPointsSource = elements.points(boundaryPointsSourceIdx,:);

% compute the normals
boundaryPointsSourceNormals = 1./sqrt(sum(elements.points(boundaryPointsSourceIdx,:).^2,2)).*elements.points(boundaryPointsSourceIdx,:); % our center is (0,0), so -> normalisation is suffices

% frequency part (robin boundary)
sourceFrequency = zeros(size(elements.points,1),1);
sourceFrequency(boundaryPointsSourceIdx) = gamma.*u1F(boundaryPointsSource(:,1),boundaryPointsSource(:,2)) + dot(squeeze(u1Fgrad(boundaryPointsSource(:,1),boundaryPointsSource(:,2))).',boundaryPointsSourceNormals.').';

% constant part (robin boundary)
sourceConstant = zeros(size(elements.points,1),1);
sourceConstant(boundaryPointsSourceIdx) = gamma.*u1C(boundaryPointsSource(:,1), boundaryPointsSource(:,2)) + dot(squeeze(u1Cgrad(boundaryPointsSource(:,1),boundaryPointsSource(:,2))).',boundaryPointsSourceNormals.').';

%%
sourceFrequency = pressure.*sourceFrequency;  
sourceConstant = pressure.*sourceConstant;
excitations = zeros(size(elements.points,1), N, 3);
excitations(:,1,1) = sourceConstant;
excitations(:,2,1) = sourceFrequency;
excitations(:,1,2) = sourceConstant;
excitations(:,2,2) = sourceFrequency;
excitations(:,1,3) = 2.*sourceConstant;
excitations(:,2,3) = 2.*sourceFrequency;
%%
[cN, U1, F] = solveWesterveltMultiLevelBoundaryExcitation(elements, omega1, beta, gamma, squeeze(kappa(:,:,1)), squeeze(excitations(:,:,1)), eta, N, N, 10^(-12));
[cN, U2, F] = solveWesterveltMultiLevelBoundaryExcitation(elements, omega2, beta, gamma, squeeze(kappa(:,:,2)), squeeze(excitations(:,:,2)), eta, N, N, 10^(-12));
[cN, U3, F] = solveWesterveltMultiLevelBoundaryExcitation(elements, omega1, beta, gamma, squeeze(kappa(:,:,3)), squeeze(excitations(:,:,3)), eta, N, N, 10^(-12));
%%
% compute the solutions on the time - space mesh

u1 = calcSolution(elements, timeMesh, squeeze(U1(N,:,:)), omega1);
u2 = calcSolution(elements, timeMesh, squeeze(U2(N,:,:)), omega2);
u3 = calcSolution(elements, timeMesh, squeeze(U3(N,:,:)), omega1);

%%
% define \Sigma our measurement manifold/discrete points
% the triangulation already defines the the edges of our doimain 1:4
% (circle sectors)

measurementEdge = 3; % positive quadrant edge

% fetch the boundary points
boundaryPointsIdx = elements.edges((elements.edges(:,3) == measurementEdge),1);

measurement_u1 = u1(:,boundaryPointsIdx);
measurement_u2 = u2(:,boundaryPointsIdx);
measurement_u3 = u3(:,boundaryPointsIdx);


