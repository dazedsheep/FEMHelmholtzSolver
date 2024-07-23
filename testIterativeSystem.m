%% Solver for the periodic Westervelt equation
clearvars

massDensity = 1000; %kg/m^3
speed_of_sound = 1480; % m/s

% signal period or center frequency
T = 10^-2;
omega = 2*pi*1/T;

% our domain
bcenter = [0,0];
brad = 1/2;
domain = [bcenter, brad];
% non linearity parameter of our domain (water = 5)
sourceValueDomain = 5;

% point scatterers and their domain
values = [0, 0];
refractionIndex = [1, 1];
% linear case
%values = [0, 0];
radii = [0.1, 0.1];
centers = [0, 0; 0.5, 0.2];

diffusivity = 10^(-9);

minHarmonics = 22; % minimum number of harmonics
nHarmonics = 22; % maximum number of harmonics

% impdeance boundary conditions --> massDensity cancels
% higher frequencies are taken into account later
beta = 1/(speed_of_sound);
gamma = 10^(-9);

meshSize = 0.005;

excitationPoints = [0;0];
% typical ultrasound pressure is 1MPa at a frequency of 1 MHz, lower
% frequency -> lower pressure!
pressure = 1*10^6;
excitationPointsSize = [0.01];
%excitationPower(1,1) = referencePressure;
%excitationPower(1,2:nHarmonics) = 0;
nExcitationHarmonics = 1;

[elements] = initializeMultiLeveLSolver(meshSize, domain);

% plot the positions of excitation(s) and source(s)
objects = getGridPointsLE(elements, [centers excitationPoints], [radii excitationPointsSize]);
figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2), objects, 'facecolor', 'interp'); shading interp;
xlabel("x [m]");
ylabel("y [m]");

% construct nonlinearity
f = constructF(elements, massDensity, speed_of_sound, refractionIndex, centers, radii, values, sourceValueDomain, true);
% construct all space dependent wave numbers for all harmonics
kappa = constructKappa(elements, diffusivity, speed_of_sound, omega, refractionIndex, centers, radii, values, nHarmonics);

% build a gaussian source
%source = 1./(speed_of_sound.^2 + 1i .* omega .* diffusivity).*referencePressure.*gaussianSource(elements, excitationPoints, 0.6);
% build a point source (regularized dirac)

% the source needs to be scaled by omega^2, this matches the SI units, we
% assume that the exication is shifted by pi
source = exp(1i.*pi).*pressure.*createPointSource(elements, excitationPoints, meshSize/2);  
excitation = zeros(size(elements.points,1),nHarmonics);
excitation(:,1) = source;

[cN, U, F] = solveWesterveltMultiLevel(elements, omega, beta, gamma, kappa, excitation, f, nHarmonics, minHarmonics, false, 10^(-12));
H = U;
U = squeeze(U(cN,:,:));
