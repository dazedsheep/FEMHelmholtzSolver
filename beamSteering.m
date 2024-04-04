clearvars

massDensity = 1000; %kg/m^3
speed_of_sound = 1480; % m/s

% signal period or center frequency
T = 10^-4;
omega = 2*pi*1/T; % angular velocity

% wavelength
lambda = speed_of_sound/(1/T);

% our domain
bcenter = [0,0];
brad = 1;
domain = [bcenter, brad];
% non linearity parameter of our domain (water = 5)
sourceValueDomain = 8;

% point scatterers and their domain
values = [0, 0];
refractionIndex = [1, 1];
% linear case
%values = [0, 0];
radii = [0.1, 0.1];
centers = [0, 0; -0.3, 0.2];

diffusivity = 10^(-9);

minHarmonics = 4; % minimum number of harmonics
nHarmonics = 4; % maximum number of harmonics

% impdeance boundary conditions --> massDensity cancels
% higher frequencies are taken into account later
beta = 1/(speed_of_sound);
gamma = 10^(-9);

meshSize = 0.005;

% build a linear array
excitationPoints = [-2*lambda/8,-lambda/8,0,lambda/8,2*lambda/8;0.8,0.8,0.8,0.8,0.8];
angles = exp(1i.*omega.*[-pi,-pi,0,pi,pi]);
% typical ultrasound pressure is 1MPa at a frequency of 1 MHz, lower
% frequency -> lower pressure!
pressure = 100*10^4;
excitationPointsSize = [0.01,0.01,0.01,0.01,0.01];
excitationPower(1,1) = pressure;
excitationPower(1,2:nHarmonics) = 0;


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

% we need to scale the reference pressure to the "point" source/linear
% array
source = pressure.*createPointSourceWithAngles(elements, excitationPoints, meshSize, angles);  
excitation = zeros(size(elements.points,1),nHarmonics);
excitation(:,1) = source;
[cN, U, F] = solveWesterveltMultiLevel(elements, omega, beta, gamma, kappa, excitation, f, nHarmonics, minHarmonics, false, 10^(-12));
H = U;
U = squeeze(U(cN,:,:));

%% array factor for linear arrays
theta = -90:0.1:90;
beamAngle = 0;
Nelements = 5;
AF = sin((Nelements.*pi.*lambda./4)/lambda.*(sin(theta) - sin(beamAngle)))./(Nelements.*((pi.*lambda./4)./lambda).*(sin(theta) - sin(beamAngle)));
