%% Solver for the periodic Westervelt equation
clearvars

massDensity = 1000; %kg/m^3
speed_of_sound = 1540; % m/s

% signal period or center frequency
T = 10^-4;
omega = 2*pi*1/T;

% our domain
bcenter = [0,0];
brad = 1;
domain = [bcenter, brad];
 
% point scatterers and their domain
values = [2.5];
refractionIndex = [1, 1];
% linear case
%values = [0, 0];
radii = [0.1];
centers = [0.0; -0.1];

diffusivity = 10^(-9);

minHarmonics = 15; % minimum number of harmonics
nHarmonics = 15; % maximum number of harmonics

% impdeance boundary conditions --> massDensity cancels
% higher frequencies are taken into account later
beta = 1/(speed_of_sound);
gamma = 1;

meshSize = 0.005;

excitationPoints = [0;0.2];
% typical ultrasound pressure is 1MPa at a frequency of 1 MHz, lower
% frequency -> lower pressure!
referencePressure = 5*10^4;
excitationPointsSize = [1];
excitationPower(1,1) = referencePressure;
excitationPower(1,2:nHarmonics) = 0;
nExcitationHarmonics = 1;

[elements] = initializeMultiLeveLSolver(meshSize, domain);

% plot the positions of excitation(s) and source(s)
%  objects = getGridPointsLE(elements, [centers excitationPoints], [radii excitationPointsSize]);
%  figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2), objects, 'facecolor', 'interp'); shading interp;
% xlabel("x [m]");
% ylabel("y [m]");


% construct nonlinearity
f = massDensity.*constructF(elements, massDensity, speed_of_sound, refractionIndex, centers, radii, values);
% construct all space dependent wave numbers for all harmonics
kappa = constructKappa(elements, diffusivity, speed_of_sound, omega, refractionIndex, centers, radii, values, nHarmonics);

% build a gaussian source
% source = -massDensity.*referencePressure.*gaussianSource(elements, excitationPoints, 0.6);
% build a point source (regularized dirac)
source = -massDensity.*speed_of_sound^2./(speed_of_sound.^2 + 1i .* omega .* diffusivity).*referencePressure.*createPointSource(elements, excitationPoints, meshSize);  
excitation = zeros(size(elements.points,1),nHarmonics);
excitation(:,1) = source;

%excitation = 1i./(speed_of_sound.^2 + 1i .* (1:nHarmonics) .* omega .* diffusivity).*excitation;
[cN, U, F] = solveWesterveltMultiLevel(elements, omega, beta, gamma, kappa, excitation, f, nHarmonics, minHarmonics, true, 10^(-12));
H = U;
U = squeeze(U(cN,:,:));
