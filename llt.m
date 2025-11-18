clearvars

% ultrasound pressure of the "point" source (in L^2(\Omega), i.e., for a single excitation frequency)

fixedError = 10^-3; % L2 error of two consecutive iterations, we stop if we reach this

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

minHarmonics = 30; % minimum number of harmonics
nHarmonics = 30; % maximum number of harmonics

gamma = 1;

excitationPoints = [0;0];
excitationSize = 0.01;

initialMeshSize = 0.0004;

speed_of_sound = 50;

% signal period or center frequency
T = 10^-2;

% create the parpool
% parpool(12);
% for this setup, the solutions start to blow up between 2000 and 3000 
pressure = 505;

omega = 2*pi* 1/T * 3;

beta = 1/speed_of_sound;
pressureIncrease = 200;
meshSize = initialMeshSize;

[elements] = initializeMultiLeveLSolver(meshSize, domain);
f = constructF(elements, massDensity, speed_of_sound, refractionIndex, centers, radii, values, sourceValueDomain, true);

kappa = constructKappa(elements, diffusivity, speed_of_sound, omega, refractionIndex, centers, radii, values, nHarmonics);

pressureSteps = 1;

for i=1:pressureSteps
    %

    % construct all space dependent wave numbers for all harmonics

    % normalise the regularised dirac
    pointSource = createPointSource(elements, excitationPoints, excitationSize);
    [area, s] = integrate_fun_trimesh(elements.opoints, elements.otri, pointSource.');
    source = -exp(1i*pi/2).*pressure.*pointSource./s;
    excitation = zeros(size(elements.points,1),nHarmonics);
    excitation(:,1) = source;
    degenerate_eta(i) = 1/max(abs(f));
    maxsource(i) = max(abs(source));

    tic
    [cN, U, G] =  solveWesterveltMultiLevelMT(elements, omega, beta, gamma, kappa, excitation, f, nHarmonics, minHarmonics, false, fixedError);
    timeMT(i) = toc;
    pressures(i) = pressure;
    for k=2:cN
        [~, int] = integrate_fun_trimesh(elements.points', elements.otri,  squeeze(abs(squeeze(real(sum(U(k,:,:)))) - squeeze(real(sum(U(k-1,:,:))))).^2 )' );
        L2errors(i,k-1) = sqrt(int);
    end
    pressure = pressure + pressureIncrease;
    numIterations(i) = cN;
    numDOF(i) = size(elements.tri,1);
    [~, int] = integrate_fun_trimesh(elements.points', elements.otri,  squeeze(abs(squeeze(real(sum(U(cN,:,:)))) - squeeze(real(sum(U(cN-1,:,:))))).^2 )' );
    error(i) = sqrt(int);
    [~, int] = integrate_fun_trimesh(elements.points', elements.otri,  squeeze(abs(squeeze(real(sum(U(cN,:,:))))).^2)' );
    L2Normt0(i) = sqrt(int);
end

%% higher speed of sound; c = 1450, B/A = 5
clearvars

% ultrasound pressure of the "point" source (in L^2(\Omega), i.e., for a single excitation frequency)

fixedError = 10^-3; % L2 error of two consecutive iterations, we stop if we reach this

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

minHarmonics = 25; % minimum number of harmonics
nHarmonics = 25; % maximum number of harmonics

gamma = 1;

excitationPoints = [0;0];
excitationSize = 0.01;

initialMeshSize = 0.0004;

speed_of_sound = 1450;

% signal period or center frequency
T = 10^-4;

% create the parpool
% parpool(12);
% for this setup, the solutions start to blow up between 2000 and 3000

beta = 1/speed_of_sound;

meshSize = initialMeshSize;

[elements] = initializeMultiLeveLSolver(meshSize, domain);

f = constructF(elements, massDensity, speed_of_sound, refractionIndex, centers, radii, values, sourceValueDomain, true);

freq = [5;7.5;10] * 10^3;
pressures = [5,7,9,11,12,13,14,15,16.2; 0.5,1.5,2.5,3.5,4.5,4.9,5.4,5.8,6.5; 0.1,0.4,0.8,1.1,1.5,1.9,2.4,2.9,3.3 ].*10^5;

for j = 1:length(freq)
    omega = 2*pi*freq(j);
    % construct all space dependent wave numbers for all harmonics
    kappa = constructKappa(elements, diffusivity, speed_of_sound, omega, refractionIndex, centers, radii, values, nHarmonics);

    for i=1:size(pressures,2)
        % normalise the regularised dirac
        pointSource = createPointSource(elements, excitationPoints, excitationSize);
        [area, s] = integrate_fun_trimesh(elements.opoints, elements.otri, pointSource.');
        source = -exp(1i*pi/2).*pressures(j,i).*pointSource./s;
        excitation = zeros(size(elements.points,1),nHarmonics);
        excitation(:,1) = source;
        % tic
        % [cN, U, F] =  solveWesterveltMultiLevel(elements, omega, beta, gamma, kappa, excitation, f, nHarmonics, minHarmonics, false, fixedError);
        % timeS = toc;
        degenerate_eta(j,i) = 1/max(abs(f));
        maxsource(j,i) = max(abs(source));
   
    end
end