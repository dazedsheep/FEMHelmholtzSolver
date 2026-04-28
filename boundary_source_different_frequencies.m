massDensity = 1000; %kg/m^3
speed_of_sound = 1480;
% signal period or center frequency of the excitation
T = 10^-5;
omega = 2*pi*1/T*sqrt(2);

% our domain
bcenter = [0,0];
brad = 0.2;
domain = [bcenter, brad];

% nonlinearity parameter of our domain (water = 5)
sourceValueDomain = 0;

% point scatterers and their domain
values = [5];
refractionIndex = [1, 1]; % this allows to adjust the speed of sound for the phantoms
% linear case
%values = [0, 0];
radii = [0.05];
% this is the diffusvity of the phantoms
diffusivityPhantoms = [20]; % this allows to adjust the diffusivity for the phantoms
centers = [0; 0];

% this is the diffusivity for the domain
diffusivity = 10^(-6);

minHarmonics = 8; % minimum number of harmonics
nHarmonics = 8; % maximum number of harmonics

% boundary paraemters

% dirichlet part
gamma = 1;

beta = 1/speed_of_sound;

meshSize = 0.0005;

% put the excitation on the boundary
excitationPoints = [0.0,0.0]; % ;-0.2,0.2 ... second source
%excitationPoints = [0.0,0.0];

% ultrasound pressure of the "point" source
pressure = 3*10^7;
excitationPointsSize = [0.001];

[elements] = initializeMultiLeveLSolver(meshSize, domain);
%%
% construct non-linearity
%f = constructF(elements, massDensity, speed_of_sound, refractionIndex, centers, radii, values, sourceValueDomain, true);
% zero nonlinearity
f = zeros(size(elements.points,1),1);

% construct all space dependent wave numbers for all harmonics
% the space dependent diffusivity is taken care of in kappa (see the
% Fourier formulation, i.e., the iteration scheme)
kappa = constructKappaS(elements, [diffusivity diffusivityPhantoms], speed_of_sound, omega, refractionIndex, centers, radii, values, nHarmonics);

% realistically piezoelectric elements are not of infinitesimal small size
source = exp(1i.*omega.*pi/2).*pressure.*createPointSourceOnBoundary(elements, excitationPoints, excitationPointsSize);  
excitation = zeros(size(elements.points,1),nHarmonics);
excitation(:,1) = source;

% solve the periodic westervelt equation with excitations on the boundary
[cN, U, F] =  solveWesterveltMultiLevelBoundaryExcitation(elements, omega, beta, gamma, kappa, excitation, f, nHarmonics, minHarmonics, false, 10^(-12));
H = U;
U = squeeze(U(cN,:,:));

P_excitation_on_the_boundary = calcPressureProfile(omega, T, H, U, cN);