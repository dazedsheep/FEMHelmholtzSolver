clearvars

% ultrasound pressure of the "point" source (in L^2(\Omega), i.e., for a single excitation frequency)

finalPressure = 10^6; % this will be the L^2 norm !!!! (\leq \delta form our Theorem)

fixedError = 10^-5; % L2 error of two consecutive iterations, we stop if we reach this

massDensity = 1000; %kg/m^3

% our domain
bcenter = [0,0];
brad = 0.05;
domain = [bcenter, brad];
% nonlinearity parameter of our domain (water = 5)
sourceValueDomain = 5;

% point scatterers and their domain
values = [0];
refractionIndex = [1];
% linear case
%values = [0, 0];
radii = [0.05];
centers = [0; 0];

diffusivity = 10^(-9);

minHarmonics = 12; % minimum number of harmonics
nHarmonics = 12; % maximum number of harmonics

gamma = 1;

excitationPoints = [0;0];
excitationSize = 0.01;

initialMeshSize = 0.0008;

speed_of_sound = 50;

% signal period or center frequency
T = 10^-2;
omega = 2*pi* 1/T;

% create the parpool
% parpool("Threads");

pressure = 1000;
frequencySteps = 1;

omega = 2*pi* 1/T;

beta = 1/speed_of_sound;

meshSize = initialMeshSize;

[elements] = initializeMultiLeveLSolver(meshSize, domain);
%
f = constructF(elements, massDensity, speed_of_sound, refractionIndex, centers, radii, values, sourceValueDomain, true);

% construct all space dependent wave numbers for all harmonics
kappa = constructKappa(elements, diffusivity, speed_of_sound, omega, refractionIndex, centers, radii, values, nHarmonics);

% normalise the regularised dirac
pointSource = createPointSource(elements, excitationPoints, excitationSize);
[area, s] = integrate_fun_trimesh(elements.opoints, elements.otri, pointSource.');
source = -exp(1i*pi/2).*pressure.*pointSource./s;
excitation = zeros(size(elements.points,1),nHarmonics);
excitation(:,1) = source;
tic
[cN, U, F] =  solveWesterveltMultiLevel(elements, omega, beta, gamma, kappa, excitation, f, nHarmonics, minHarmonics, false, fixedError);
timeS = toc;

tic
[cN, A, G] =  solveWesterveltMultiLevelMT(elements, omega, beta, gamma, kappa, excitation, f, nHarmonics, minHarmonics, false, fixedError);
timeMT = toc;

