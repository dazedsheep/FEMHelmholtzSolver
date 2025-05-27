%% Try to find some thresholds for $\delta$ from our Theorem on the existence of solutions to the nonlinear Westervelt equation
% clearly $\delta$ depends on the speed of sound, frequency, domain

%% c = 50, f = 100 Hz, B/A = 5, fix the error and start with a coarse mesh and gradually increase the pressure
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

initialMeshSize = 0.05;

speed_of_sound = 1450;

% signal period or center frequency
T = 10^-4;

% create the parpool
pressure = 1*10^5;

omega = 2*pi* 1/T;

pressureIncrease = 10000;

pressureSteps = 20;
beta = 1/speed_of_sound;

for i = 1:pressureSteps


    meshSize(i) = initialMeshSize;

    while true
        [elements] = initializeMultiLeveLSolver(meshSize(i), domain);


    meshSize(i) = initialMeshSize;

    while true
        [elements] = initializeMultiLeveLSolver(meshSize(i), domain);

        f = constructF(elements, massDensity, speed_of_sound, refractionIndex, centers, radii, values, sourceValueDomain, true);
        f = constructF(elements, massDensity, speed_of_sound, refractionIndex, centers, radii, values, sourceValueDomain, true);

        % construct all space dependent wave numbers for all harmonics
        kappa = constructKappa(elements, diffusivity, speed_of_sound, omega, refractionIndex, centers, radii, values, nHarmonics);
        % construct all space dependent wave numbers for all harmonics
        kappa = constructKappa(elements, diffusivity, speed_of_sound, omega, refractionIndex, centers, radii, values, nHarmonics);

        % normalise the regularised dirac
        pointSource = createPointSource(elements, excitationPoints, excitationSize);
        [area, s] = integrate_fun_trimesh(elements.opoints, elements.otri, pointSource.');
        source = -exp(1i*pi/2).*pressure.*pointSource./s;
        excitation = zeros(size(elements.points,1),nHarmonics);
        excitation(:,1) = source;
        degenerate_eta(i) = 1/max(abs(f));
        maxsource(i) = max(abs(source));
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
        numIterations(i) = cN;
        numDOF(i) = size(elements.tri,1);
        [~, int] = integrate_fun_trimesh(elements.points', elements.otri,  squeeze(abs(squeeze(real(sum(U(cN,:,:)))) - squeeze(real(sum(U(cN-1,:,:))))).^2 )' );
        error(i) = sqrt(int);
        [~, int] = integrate_fun_trimesh(elements.points', elements.otri,  squeeze(abs(squeeze(real(sum(U(cN,:,:))))).^2)' );
        L2Normt0(i) = sqrt(int);
        if L2errors(i,cN-1) < fixedError
            break;
        end
        % if we did not reach the error level increase the resolution
        meshSize(i) = meshSize(i) / 2;
    end
        tic
        [cN, U, G] =  solveWesterveltMultiLevelMT(elements, omega, beta, gamma, kappa, excitation, f, nHarmonics, minHarmonics, false, fixedError);
        timeMT(i) = toc;
        pressures(i) = pressure;
        for k=2:cN
            [~, int] = integrate_fun_trimesh(elements.points', elements.otri,  squeeze(abs(squeeze(real(sum(U(k,:,:)))) - squeeze(real(sum(U(k-1,:,:))))).^2 )' );
            L2errors(i,k-1) = sqrt(int);
        end
        numIterations(i) = cN;
        numDOF(i) = size(elements.tri,1);
        [~, int] = integrate_fun_trimesh(elements.points', elements.otri,  squeeze(abs(squeeze(real(sum(U(cN,:,:)))) - squeeze(real(sum(U(cN-1,:,:))))).^2 )' );
        error(i) = sqrt(int);
        [~, int] = integrate_fun_trimesh(elements.points', elements.otri,  squeeze(abs(squeeze(real(sum(U(cN,:,:))))).^2)' );
        L2Normt0(i) = sqrt(int);
        if L2errors(i,cN-1) < fixedError
            break;
        end
        % if we did not reach the error level increase the resolution
        meshSize(i) = meshSize(i) / 2;
    end
    i
    pressure = pressure + pressureIncrease;


end

%% c = 50, f = 100 Hz, B/A = 5
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

initialMeshSize = 0.0005;

speed_of_sound = 50;

% signal period or center frequency
T = 10^-2;

% create the parpool
% parpool(12);
% for this setup, the solutions start to blow up between 2000 and 3000
% for this setup, the solutions start to blow up between 2000 and 3000
pressure = 3000;

omega = 2*pi* 1/T;

beta = 1/speed_of_sound;
pressureIncrease = 200;
meshSize = initialMeshSize;

[elements] = initializeMultiLeveLSolver(meshSize, domain);

pressureSteps = 7;

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
    degenerate_eta(i) = 1/max(abs(f));
    maxsource(i) = max(abs(source));

    tic
    [cN, U, G] =  solveWesterveltMultiLevelMT(elements, omega, beta, gamma, kappa, excitation, f, nHarmonics, minHarmonics, false, 0);
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
%%
%%
figure,
hold on
for i = 1:pressureSteps
    plot(log10(L2errors(i,:)))
end
hold off

sc = 7;

figure, scatter(pressures(1:sc), L2Normt0(1:sc))
hold on
fit = (L2Normt0(2) - L2Normt0(1))/(pressures(2)-pressures(1)) * (0:(pressureSteps-1))*pressureIncrease + L2Normt0(1);
plot(pressures(1:sc),fit(1:sc));
hold off
%% c = 50, f = 100 Hz, B/A = 5
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

initialMeshSize = 0.0005;

speed_of_sound = 50;

% signal period or center frequency
T = 10^-2;

% create the parpool
% parpool(12);
% for this setup, the solutions start to blow up between 2000 and 3000
% for this setup, the solutions start to blow up between 2000 and 3000
pressure = 1000;

omega = 2*pi* 1/T;

beta = 1/speed_of_sound;
pressureIncrease = 200;
meshSize = initialMeshSize;

[elements] = initializeMultiLeveLSolver(meshSize, domain);

pressureSteps = 17;

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
%%
%%
figure,
hold on
for i = 1:pressureSteps
    plot(log10(L2errors(i,:)))
end
hold off

sc = 17;

figure, scatter(pressures(1:sc), L2Normt0(1:sc))
hold on
fit = (L2Normt0(2) - L2Normt0(1))/(pressures(2)-pressures(1)) * (0:(pressureSteps-1))*pressureIncrease + L2Normt0(1);
plot(pressures(1:sc),fit(1:sc));
hold off

%% c = 50, f = 200 Hz B/A = 5
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

initialMeshSize = 0.0005;

speed_of_sound = 50;

% signal period or center frequency
T = 10^-2;

% create the parpool
% parpool(12);
% for this setup, the solutions start to blow up between 2000 and 3000
pressure = 500;

omega = 2*pi* 1/T * 2;

beta = 1/speed_of_sound;
pressureIncrease = 100;
meshSize = initialMeshSize;

[elements] = initializeMultiLeveLSolver(meshSize, domain);

pressureSteps = 7;

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
    degenerate_eta(i) = 1/max(abs(f));
    maxsource(i) = max(abs(source));

    tic
    [cN, U, G] =  solveWesterveltMultiLevelMT(elements, omega, beta, gamma, kappa, excitation, f, nHarmonics, minHarmonics, false, 0);
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
%%
%%
figure,
hold on
for i = 1:pressureSteps
    plot(log10(L2errors(i,:)))
end
hold off

sc = 7;

figure, scatter(pressures(1:sc), L2Normt0(1:sc))
hold on
fit = (L2Normt0(2) - L2Normt0(1))/(pressures(2)-pressures(1)) * (0:(pressureSteps-1))*pressureIncrease + L2Normt0(1);
plot(pressures(1:sc),fit(1:sc));
hold off
%% c = 50, f = 300 Hz
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

initialMeshSize = 0.0005;

speed_of_sound = 50;

% signal period or center frequency
T = 10^-2;

% create the parpool
% parpool(12);
% for this setup, the solutions start to blow up between 2000 and 3000
pressure = 200;

omega = 2*pi* 1/T * 3;

beta = 1/speed_of_sound;
pressureIncrease = 30;
meshSize = initialMeshSize;

[elements] = initializeMultiLeveLSolver(meshSize, domain);

pressureSteps = 7;

for i=1:pressureSteps
    %
    f = constructF(elements, massDensity, speed_of_sound, refractionIndex, centers, radii, values, sourceValueDomain, true);

    % construct all space dependent wave numbers for all harmonics
    kappa = constructKappa(elements, diffusivity, speed_of_sound, omega, refractionIndex, centers, radii, values, nHarmonics);

    % normalise the regularised dirac
    pointSource = createPointSource(elements, excitationPoints, excitationSize);
    [area, s] = integrate_fun_trimesh(elements.opoints, elements.otri, pointSource.');
    [area, s] = integrate_fun_trimesh(elements.opoints, elements.otri, pointSource.');
    source = -exp(1i*pi/2).*pressure.*pointSource./s;
    degenerate_eta(i) = 1/max(abs(f));
    maxsource(i) = max(abs(source));
    excitation = zeros(size(elements.points,1),nHarmonics);
    excitation(:,1) = source;
    % tic
    % [cN, U, F] =  solveWesterveltMultiLevel(elements, omega, beta, gamma, kappa, excitation, f, nHarmonics, minHarmonics, false, fixedError);
    % timeS = toc;

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
%%
%%
figure,
hold on
for i = 1:pressureSteps
    plot(log10(L2errors(i,:)))
end
hold off

sc = 7;

figure, scatter(pressures(1:sc), L2Normt0(1:sc))
hold on
fit = (L2Normt0(2) - L2Normt0(1))/(pressures(2)-pressures(1)) * (0:(pressureSteps-1))*pressureIncrease + L2Normt0(1);
plot(pressures(1:sc),fit(1:sc));
hold off

%% higher speed of sound; c = 1450, f = 5 kHz, B/A = 5
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
pressure = 5*10^5;

omega = 2*pi* 1/T * 1/2;

beta = 1/speed_of_sound;
pressureIncrease = 200000;
meshSize = initialMeshSize;

[elements] = initializeMultiLeveLSolver(meshSize, domain);

f = constructF(elements, massDensity, speed_of_sound, refractionIndex, centers, radii, values, sourceValueDomain, true);

% construct all space dependent wave numbers for all harmonics
kappa = constructKappa(elements, diffusivity, speed_of_sound, omega, refractionIndex, centers, radii, values, nHarmonics);
pressures = [5,7,9,11,12,14,15,16.2].*10^5;


for i=1:length(pressures)
    % normalise the regularised dirac
    pointSource = createPointSource(elements, excitationPoints, excitationSize);
    [area, s] = integrate_fun_trimesh(elements.opoints, elements.otri, pointSource.');
    source = -exp(1i*pi/2).*pressures(i).*pointSource./s;
    excitation = zeros(size(elements.points,1),nHarmonics);
    excitation(:,1) = source;
    % tic
    % [cN, U, F] =  solveWesterveltMultiLevel(elements, omega, beta, gamma, kappa, excitation, f, nHarmonics, minHarmonics, false, fixedError);
    % timeS = toc;
    degenerate_eta(i) = 1/max(abs(f));
    maxsource(i) = max(abs(source));
    tic
    [cN, U, G] =  solveWesterveltMultiLevelMT(elements, omega, beta, gamma, kappa, excitation, f, nHarmonics, minHarmonics, false, fixedError);
    timeMT(i) = toc;

    for k=2:cN
        [~, int] = integrate_fun_trimesh(elements.points', elements.otri,  squeeze(abs(squeeze(real(sum(U(k,:,:)))) - squeeze(real(sum(U(k-1,:,:))))).^2 )' );
        L2errors(i,k-1) = sqrt(int);
    end
    numIterations(i) = cN;
    numDOF(i) = size(elements.tri,1);
    [~, int] = integrate_fun_trimesh(elements.points', elements.otri,  squeeze(abs(squeeze(real(sum(U(cN,:,:)))) - squeeze(real(sum(U(cN-1,:,:))))).^2 )' );
    error(i) = sqrt(int);
    [~, int] = integrate_fun_trimesh(elements.points', elements.otri,  squeeze(abs(squeeze(real(sum(U(cN,:,:))))).^2)' );
    L2Normt0(i) = sqrt(int);
    degen(i)= max(max(abs(f.* squeeze(abs(squeeze(real(sum(U(cN,:,:)))))))));
end
%%
figure,
hold on
for i = 1:length(pressures)
    plot(log10(L2errors(i,:)))
end
hold off

sc = length(pressures);

figure, scatter(pressures(1:sc), L2Normt0(1:sc))
hold on
fit = (L2Normt0(2) - L2Normt0(1))/(pressures(2)-pressures(1)) * (pressures(1:sc) - pressures(1)) + L2Normt0(1);
plot(pressures(1:sc),fit(1:sc));
hold off
%% higher speed of sound; c = 1450, f = 7.5 kHz, B/A = 5
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

% create the parpool
% parpool(12);
% for this setup, the solutions start to blow up between 2000 and 3000

omega = 2*pi* 7.5 * 10^3;

beta = 1/speed_of_sound;
meshSize = initialMeshSize;

[elements] = initializeMultiLeveLSolver(meshSize, domain);

f = constructF(elements, massDensity, speed_of_sound, refractionIndex, centers, radii, values, sourceValueDomain, true);

% construct all space dependent wave numbers for all harmonics
kappa = constructKappa(elements, diffusivity, speed_of_sound, omega, refractionIndex, centers, radii, values, nHarmonics);
pressures = [0.5,1.5,2.5,3.5,4.5,4.9,5.4,5.8,6.5].*10^5;

for i=1:length(pressures)
    % normalise the regularised dirac
    pointSource = createPointSource(elements, excitationPoints, excitationSize);
    [area, s] = integrate_fun_trimesh(elements.opoints, elements.otri, pointSource.');
    source = -exp(1i*pi/2).*pressures(i).*pointSource./s;
    excitation = zeros(size(elements.points,1),nHarmonics);
    excitation(:,1) = source;
    % tic
    % [cN, U, F] =  solveWesterveltMultiLevel(elements, omega, beta, gamma, kappa, excitation, f, nHarmonics, minHarmonics, false, fixedError);
    % timeS = toc;
    degenerate_eta(i) = 1/max(abs(f));
    maxsource(i) = max(abs(source));
    tic
    [cN, U, G] =  solveWesterveltMultiLevelMT(elements, omega, beta, gamma, kappa, excitation, f, nHarmonics, minHarmonics, false, fixedError);
    timeMT(i) = toc;

    for k=2:cN
        [~, int] = integrate_fun_trimesh(elements.points', elements.otri,  squeeze(abs(squeeze(real(sum(U(k,:,:)))) - squeeze(real(sum(U(k-1,:,:))))).^2 )' );
        L2errors(i,k-1) = sqrt(int);
    end
    numIterations(i) = cN;
    numDOF(i) = size(elements.tri,1);
    [~, int] = integrate_fun_trimesh(elements.points', elements.otri,  squeeze(abs(squeeze(real(sum(U(cN,:,:)))) - squeeze(real(sum(U(cN-1,:,:))))).^2 )' );
    error(i) = sqrt(int);
    [~, int] = integrate_fun_trimesh(elements.points', elements.otri,  squeeze(abs(squeeze(real(sum(U(cN,:,:))))).^2)' );
    L2Normt0(i) = sqrt(int);
    degen(i)= max(max(abs(f.* squeeze(abs(squeeze(real(sum(U(cN,:,:)))))))));
end
%%
figure,
hold on
for i = 1:length(pressures)
    plot(log10(L2errors(i,:)))
end
hold off

sc = length(pressures);

figure, scatter(pressures(1:sc), L2Normt0(1:sc))
hold on
fit = (L2Normt0(2) - L2Normt0(1))/(pressures(2)-pressures(1)) * (pressures(1:sc) - pressures(1)) + L2Normt0(1);
plot(pressures(1:sc),fit(1:sc));
hold off
%%
%% higher speed of sound; c = 1450, f = 10 kHz, B/A = 5
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

omega = 2*pi* 10^4;

beta = 1/speed_of_sound;
meshSize = initialMeshSize;

[elements] = initializeMultiLeveLSolver(meshSize, domain);

f = constructF(elements, massDensity, speed_of_sound, refractionIndex, centers, radii, values, sourceValueDomain, true);

% construct all space dependent wave numbers for all harmonics
kappa = constructKappa(elements, diffusivity, speed_of_sound, omega, refractionIndex, centers, radii, values, nHarmonics);
pressures = [0.4,0.8,1.1,1.5,1.9,2.4,2.9,3.4].*10^5;

for i=1:length(pressures)
    % normalise the regularised dirac
    pointSource = createPointSource(elements, excitationPoints, excitationSize);
    [area, s] = integrate_fun_trimesh(elements.opoints, elements.otri, pointSource.');
    source = -exp(1i*pi/2).*pressures(i).*pointSource./s;
    excitation = zeros(size(elements.points,1),nHarmonics);
    excitation(:,1) = source;
    % tic
    % [cN, U, F] =  solveWesterveltMultiLevel(elements, omega, beta, gamma, kappa, excitation, f, nHarmonics, minHarmonics, false, fixedError);
    % timeS = toc;
    degenerate_eta(i) = 1/max(abs(f));
    maxsource(i) = max(abs(source));
    tic
    [cN, U, G] =  solveWesterveltMultiLevelMT(elements, omega, beta, gamma, kappa, excitation, f, nHarmonics, minHarmonics, false, fixedError);
    timeMT(i) = toc;

    for k=2:cN
        [~, int] = integrate_fun_trimesh(elements.points', elements.otri,  squeeze(abs(squeeze(real(sum(U(k,:,:)))) - squeeze(real(sum(U(k-1,:,:))))).^2 )' );
        L2errors(i,k-1) = sqrt(int);
    end
    numIterations(i) = cN;
    numDOF(i) = size(elements.tri,1);
    [~, int] = integrate_fun_trimesh(elements.points', elements.otri,  squeeze(abs(squeeze(real(sum(U(cN,:,:)))) - squeeze(real(sum(U(cN-1,:,:))))).^2 )' );
    error(i) = sqrt(int);
    [~, int] = integrate_fun_trimesh(elements.points', elements.otri,  squeeze(abs(squeeze(real(sum(U(cN,:,:))))).^2)' );
    L2Normt0(i) = sqrt(int);
    degen(i)= max(max(abs(f.* squeeze(abs(squeeze(real(sum(U(cN,:,:)))))))));
end
%%
figure,
hold on
for i = 1:length(pressures)
    plot(log10(L2errors(i,:)))
end
hold off

% figure,
% hold on
% for i = 1:pressureSteps
% 
%     plot(abs(L2errors(i,:) - L2errors(i,:)))
% end
% hold off

% figure,
% hold on
% for i = 1:pressureSteps
% 
%     plot(abs(L2errors(i,:) - L2errors(i,:)))
% end
% hold off

sc = length(pressures);

figure, scatter(pressures(1:sc), L2Normt0(1:sc))
hold on
fit = (L2Normt0(2) - L2Normt0(1))/(pressures(2)-pressures(1)) * (pressures(1:sc) - pressures(1)) + L2Normt0(1);
plot(pressures(1:sc),fit(1:sc));
hold off

%%
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
        tic
        [cN, U, G] =  solveWesterveltMultiLevelMT(elements, omega, beta, gamma, kappa, excitation, f, nHarmonics, minHarmonics, false, fixedError);
        timeMT(j,i) = toc;

        for k=2:cN
            [~, int] = integrate_fun_trimesh(elements.points', elements.otri,  squeeze(abs(squeeze(real(sum(U(k,:,:)))) - squeeze(real(sum(U(k-1,:,:))))).^2 )' );
            L2errors(j,i,k-1) = sqrt(int);
        end
        numIterations(j,i) = cN;
        numDOF(j,i) = size(elements.tri,1);
        [~, int] = integrate_fun_trimesh(elements.points', elements.otri,  squeeze(abs(squeeze(real(sum(U(cN,:,:)))) - squeeze(real(sum(U(cN-1,:,:))))).^2 )' );
        error(j,i) = sqrt(int);
        [~, int] = integrate_fun_trimesh(elements.points', elements.otri,  squeeze(abs(squeeze(real(sum(U(cN,:,:))))).^2)' );
        L2Normt0(j,i) = sqrt(int);
        degen(j,i)= max(max(abs(f.* squeeze(abs(squeeze(real(sum(U(cN,:,:)))))))));
    end
end
%% plot all the results
sc = size(pressures,2);
legindex = 1;
figure, 
for j = 1:length(freq)
    scatter(pressures(j,1:sc), L2Normt0(j,1:sc))
    leg{legindex} = sprintf('f = %d kHz', freq(j)./10^3);
    hold on
    fitlin = (L2Normt0(j,2) - L2Normt0(j,1))/(pressures(j,2)-pressures(j,1)) * (pressures(j,1:sc) - pressures(j,1)) + L2Normt0(j,1);
    plot(pressures(j,1:sc),fitlin(1:sc));
    leg{legindex + 1} = sprintf('linear fit (before saturation), f = %d kHz', freq(j)./10^3);
    leg{legindex + 2} = sprintf('poly fit (degree = 3), f = %d kHz', freq(j)./10^3);

    legindex = legindex + 3;
    fitnl = polyfit(pressures(j,1:(sc-1)), L2Normt0(j,1:(sc-1)),3);
    plot(pressures(j,1:(sc-1)),polyval(fitnl,pressures(j,1:(sc-1))));
    hold on
end
xlabel('$||\hat{h}||_{L^2(\Omega)}$', 'Interpreter','latex');
ylabel('$||p^N(0,\cdot)||_{L^2(\Omega)}$', 'Interpreter', 'latex');
legend(leg);
hold off
%%
for j = 1:length(freq)
    figure,
    hold on
    for i = 1:length(pressures)
        plot(2:25,log10(squeeze(L2errors(j,i,:))))
        errorleg{i} = sprintf('$||\\hat{h}||_{L^2(\\Omega)}$ = %d', pressures(j,i));
    end
    ylabel('$L^2$error log10$(e_j(0)))$', 'Interpreter','latex');
    xlabel('Iterations (N)');
    ylim([-8,8]);
    legend(errorleg, 'Interpreter', 'latex');
    hold off
end

%%
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
        tic
        [cN, U, G] =  solveWesterveltMultiLevelMT(elements, omega, beta, gamma, kappa, excitation, f, nHarmonics, minHarmonics, false, fixedError);
        timeMT(j,i) = toc;

        for k=2:cN
            [~, int] = integrate_fun_trimesh(elements.points', elements.otri,  squeeze(abs(squeeze(real(sum(U(k,:,:)))) - squeeze(real(sum(U(k-1,:,:))))).^2 )' );
            L2errors(j,i,k-1) = sqrt(int);
        end
        numIterations(j,i) = cN;
        numDOF(j,i) = size(elements.tri,1);
        [~, int] = integrate_fun_trimesh(elements.points', elements.otri,  squeeze(abs(squeeze(real(sum(U(cN,:,:)))) - squeeze(real(sum(U(cN-1,:,:))))).^2 )' );
        error(j,i) = sqrt(int);
        [~, int] = integrate_fun_trimesh(elements.points', elements.otri,  squeeze(abs(squeeze(real(sum(U(cN,:,:))))).^2)' );
        L2Normt0(j,i) = sqrt(int);
        degen(j,i)= max(max(abs(f.* squeeze(abs(squeeze(real(sum(U(cN,:,:)))))))));
    end
end
%% plot all the results
sc = size(pressures,2);
legindex = 1;
figure, 
for j = 1:length(freq)
    scatter(pressures(j,1:sc), L2Normt0(j,1:sc))
    leg{legindex} = sprintf('f = %d kHz', freq(j)./10^3);
    hold on
    fitlin = (L2Normt0(j,2) - L2Normt0(j,1))/(pressures(j,2)-pressures(j,1)) * (pressures(j,1:sc) - pressures(j,1)) + L2Normt0(j,1);
    plot(pressures(j,1:sc),fitlin(1:sc));
    leg{legindex + 1} = sprintf('linear fit (before saturation), f = %d kHz', freq(j)./10^3);
    leg{legindex + 2} = sprintf('poly fit (degree = 3), f = %d kHz', freq(j)./10^3);

    legindex = legindex + 3;
    fitnl = polyfit(pressures(j,1:(sc-1)), L2Normt0(j,1:(sc-1)),3);
    plot(pressures(j,1:(sc-1)),polyval(fitnl,pressures(j,1:(sc-1))));
    hold on
end
xlabel('$||\hat{h}||_{L^2(\Omega)}$', 'Interpreter','latex');
ylabel('$||p^N(0,\cdot)||_{L^2(\Omega)}$', 'Interpreter', 'latex');
legend(leg);
hold off
%%
for j = 1:length(freq)
    figure,
    hold on
    for i = 1:length(pressures)
        plot(2:25,log10(squeeze(L2errors(j,i,:))))
        errorleg{i} = sprintf('$||\\hat{h}||_{L^2(\\Omega)}$ = %d', pressures(j,i));
    end
    ylabel('$L^2$error log10$(e_j(0)))$', 'Interpreter','latex');
    xlabel('Iterations (N)');
    ylim([-8,8]);
    legend(errorleg, 'Interpreter', 'latex');
    hold off
end

%% do all the simulations at once for B/A=5,
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

initialMeshSize = 0.0005;

speed_of_sound = 50;

% signal period or center frequency
T = 10^-2;

frequencySteps = 1;

beta = 1/speed_of_sound;

meshSize = initialMeshSize;

[elements] = initializeMultiLeveLSolver(meshSize, domain);

pressureSteps=7;
startPressure = [3000,800,325];
pressureInc = [200, 100, 35];
for j = 1:3 % 100 Hz, 200 Hz, 300 Hz
    pressure = startPressure(j);
    omega = 2*pi* 1/T * j;
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
%        tic
        [cN, U, G] =  solveWesterveltMultiLevelMT(elements, omega, beta, gamma, kappa, excitation, f, nHarmonics, minHarmonics, false, 0);
%        timeMT(i) = toc;
        pressures(j,i) = pressure;
        for k=2:cN
            [~, int] = integrate_fun_trimesh(elements.points', elements.otri,  squeeze(abs(squeeze(real(sum(U(k,:,:)))) - squeeze(real(sum(U(k-1,:,:))))).^2 )' );
            L2errors(j,i,k-1) = sqrt(int);
        end
        pressure = pressure + pressureInc(j);
        numIterations(j,i) = cN;
        numDOF(j,i) = size(elements.tri,1);
        [~, int] = integrate_fun_trimesh(elements.points', elements.otri,  squeeze(abs(squeeze(real(sum(U(cN,:,:)))) - squeeze(real(sum(U(cN-1,:,:))))).^2 )' );
        error(j,i) = sqrt(int);
        [~, int] = integrate_fun_trimesh(elements.points', elements.otri,  squeeze(abs(squeeze(real(sum(U(cN,:,:))))).^2)' );
        L2Normt0(j,i) = sqrt(int);
        % compute the degeneration condition
        degenCond(j,i) = max(max(abs(f.*real(sum(U(cN,:,:))))));
        
    end
end
%% plot things
for j=1:3
    figure,
    hold on
    for i = 1:pressureSteps
        plot(log10(squeeze(L2errors(j,i,:))))
        leg = sprintf('$||\\hat{h}||_{L^2(\\Omega)}$ = %d', pressures(j,i));
        mylegend{i} = leg;
    end
    ylabel('$L^2$ error $log10(e_j(0))$', 'Interpreter', 'latex');
    xlabel('Iterations (N)', 'Interpreter', 'latex');
    legend(mylegend, 'Interpreter','latex');
    hold off
end
%% and plot the energy for different frequencies
sc = 7;

figure,
for j=1:3
    scatter(pressures(j,1:sc), L2Normt0(j,1:sc))
    hold on
    fit = (L2Normt0(j,2) - L2Normt0(j,1))/(pressures(j,2)-pressures(j,1)) * (0:(pressureSteps-1))*pressureInc(j) + L2Normt0(j,1);
    plot(pressures(j,1:sc),fit(1:sc));
    hold on
end
ylabel('$||p^N(0,\cdot)||_{L^2(\Omega)}$','Interpreter','latex');
xlabel('$||\hat{h}||_{L^2(\Omega)}$', 'Interpreter','latex');
legend('f = 100 Hz', '$ c_1 ||\hat{h}||_{L^2(\Omega)}, f = 100 Hz$', 'f = 200 Hz', '$c_2 ||\hat{h}||_{L^2(\Omega)}, f = 200 Hz$', 'f = 300 Hz','$c_3 ||\hat{h}||_{L^2(\Omega)}$, f = 300 Hz', 'Interpreter','latex');
%%
[area, int] = integrate_fun_trimesh(elements.points', elements.otri,  squeeze(abs(squeeze(real(sum(U(cN,:,:))))).^2)' );
c1 = L2Normt0(1,1)/pressures(1,1);
fit = (L2Normt0(1,2) - L2Normt0(1,1))/(pressures(1,2)-pressures(1,1)) * (0:(pressureSteps-1))*pressureInc(1) + L2Normt0(1,1);
1/(c1 * 2* 1/100 * sqrt(area)* max(abs(f)))

%% 120 kPa, c = 150 m/s, f = 100 Hz;
massdensity = 1000;
ba = 5;
c = [50, 50, 50; 75, 75, 75; 100, 100, 100];
x = [4000, 1400, 500; 18000, 6000, 3000 ;40000, 18000, 9000];
f = [100, 200, 300; 100, 200, 300 ;100, 200, 300];
eta = (1+ 1./2.*ba)./(massDensity.*(c));

omegas = 2*pi*f;
kappas = omegas./c;
% plot kappas, eta, x
