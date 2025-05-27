%% find numerical values for delta using a greedy algorithm
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

minHarmonics = 5; % minimum number of harmonics
nHarmonics = 30; % maximum number of harmonics

gamma = 1;

excitationPoints = [0;0];
excitationSize = 0.01;

initialMeshSize = 0.005;


% create the parpool
% parpool(12);
% for this setup, the solutions start to blow up between 2000 and 3000 



pressureIncrease = 10000;
meshSize = initialMeshSize;

speed_of_sound = 300:100:1500;
freq = 1000:500:(10^5);
pressure = 400:400:10^6;

[C,F,P] = meshgrid(speed_of_sound, freq, pressure);

[elements] = initializeMultiLeveLSolver(meshSize, domain);
ind= zeros(size(C));
for i = 1:size(C,2)
    for j = 1: size(F,1)
        omega = 2*pi*freq(j);
        kappa = constructKappa(elements, diffusivity, speed_of_sound(i), omega, refractionIndex, centers, radii, values, nHarmonics);
        if (nHarmonics*meshSize * abs(kappa(1,1)) ) < 1
            ind(i,j,:) = 1;
        end
    end
end

%%
for k = 1:size(C,2)
    for j = 1:size(F,1)
        if (ind(k,j,1) == 0)
            break;
        end
        for i = 1:size(P,3)
            
        end
    end
end

%%
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
    L2Normt0(i) = sqrt(int);thi
end

%% 120 kPa, c = 150 m/s, f = 100 Hz;
massdensity = 1000;
ba = 5;
c = [50, 50, 50, 100, 100, 100];
x = [4000, 1400, 500, 40000, 18000, 9000];
f = [100, 200, 300, 100, 200, 300];
eta = (1+ 1./2.*ba)./(massDensity.*(c));

omegas = 2*pi*f;
kappas = omegas./c;
% plot kappas, eta, x