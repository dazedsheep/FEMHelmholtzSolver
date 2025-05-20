%% Try to find some thresholds for $\delta$ from our Theorem on the existence of solutions to the nonlinear Westervelt equation
% clearly $\delta$ depends on the speed of sound, frequency, domain
clearvars

% ultrasound pressure of the "point" source (in L^2(\Omega), i.e., for a single excitation frequency)

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
nHarmonics = 20; % maximum number of harmonics

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

pressure = 10;
frequencySteps = 1;

omega = 2*pi* 1/T;

beta = 1/speed_of_sound;

meshSize = initialMeshSize;

[elements] = initializeMultiLeveLSolver(meshSize, domain);

pressureSteps=5;

for i=1:pressureSteps
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
% tic
% [cN, U, F] =  solveWesterveltMultiLevel(elements, omega, beta, gamma, kappa, excitation, f, nHarmonics, minHarmonics, false, fixedError);
% timeS = toc;

tic
[cN, A, G] =  solveWesterveltMultiLevelMT(elements, omega, beta, gamma, kappa, excitation, f, nHarmonics, minHarmonics, false, fixedError);
timeMT(i) = toc;
pressures(i) = pressure;
pressure = pressure * 5;

numIterations(i) = cN;
numDOF(i) = size(elements.tri,1);
[~, int] = integrate_fun_trimesh(elements.points', elements.otri,  squeeze(abs(squeeze(real(sum(A(cN,:,:)))) - squeeze(real(sum(A(cN-1,:,:))))).^2 )' );
error(i) = sqrt(int);

end

%%
clearvars

% ultrasound pressure of the "point" source (in L^2(\Omega), i.e., for a single excitation frequency)

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
nHarmonics = 30; % maximum number of harmonics

gamma = 1;

excitationPoints = [0;0];
excitationSize = 0.01;

initialMeshSize = 0.0008;

speed_of_sound = 50;

% signal period or center frequency
T = 10^-2;

% create the parpool
% parpool(12);
% for this setup, the solutions start to blow up between 2000 and 3000 
pressure = 2000;
frequencySteps = 1;

omega = 2*pi* 1/T;

beta = 1/speed_of_sound;

meshSize = initialMeshSize;

[elements] = initializeMultiLeveLSolver(meshSize, domain);

pressureSteps=15;

for i=1:pressureSteps
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
% tic
% [cN, U, F] =  solveWesterveltMultiLevel(elements, omega, beta, gamma, kappa, excitation, f, nHarmonics, minHarmonics, false, fixedError);
% timeS = toc;

tic
[cN, U, G] =  solveWesterveltMultiLevelMT(elements, omega, beta, gamma, kappa, excitation, f, nHarmonics, minHarmonics, false, fixedError);
timeMT(i) = toc;
pressures(i) = pressure;
for k=2:minHarmonics
        [~, int] = integrate_fun_trimesh(elements.points', elements.otri,  squeeze(abs(squeeze(real(sum(U(k,:,:)))) - squeeze(real(sum(U(k-1,:,:))))).^2 )' );
        L2errors(i,k-1) = sqrt(int);    
    end
pressure = pressure + 100;
numIterations(i) = cN;
numDOF(i) = size(elements.tri,1);
[~, int] = integrate_fun_trimesh(elements.points', elements.otri,  squeeze(abs(squeeze(real(sum(U(cN,:,:)))) - squeeze(real(sum(U(cN-1,:,:))))).^2 )' );
error(i) = sqrt(int);

end