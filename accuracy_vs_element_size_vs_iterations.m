%% runtime evals

clearvars

massDensity = 1000; %kg/m^3
speed_of_sound = 1450; % m/s

% signal period or center frequency 10 kHz
T = 10^-5;
omega = 2*pi*1/T;

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

minHarmonics = 20; % minimum number of harmonics
nHarmonics = 20; % maximum number of harmonics

gamma = 1;

beta = 1/speed_of_sound;

meshSteps = 6;

initialMeshsize = 0.0025;

meshdivfactor = 2;
x = (-0.05):(initialMeshsize/((meshdivfactor)^(meshSteps+1))):(0.05);
y = x;

%
excitationPoints = [0;0];  
%
excitationSize = 0.01;
pressure = 2.5*10^4;
% figure, trisurf(source_elements.tri(:,1:3), source_elements.points(:,1), source_elements.points(:,2), source, 'facecolor', 'interp'); shading interp;
%%
for i = 1:meshSteps

    meshSize = initialMeshsize/(meshdivfactor^(i-1));
    [elements] = initializeMultiLeveLSolver(meshSize, domain);
    f = constructF(elements, massDensity, speed_of_sound, refractionIndex, centers, radii, values, sourceValueDomain, true);
    % construct all space dependent wave numbers for all harmonics
    kappa = constructKappa(elements, diffusivity, speed_of_sound, omega, refractionIndex, centers, radii, values, nHarmonics);
    pointSource = createPointSource(elements, excitationPoints, excitationSize);
    [area, s] = integrate_fun_trimesh(elements.opoints, elements.otri, pointSource.');
    source = -exp(1i*pi/2).*pressure.*pointSource./s;  
    excitation = zeros(size(elements.points,1),nHarmonics);
    excitation(:,1) = source;
    tic
    [cN, U, F] =  solveWesterveltMultiLevelMT(elements, omega, beta, gamma, kappa, excitation, f, nHarmonics, minHarmonics, false, 0);
    solveTime(i) = toc;
    numIterations(i) = cN;
    H = U;
    % calculate the difference of consecutive solution over the iterations
    for k=2:minHarmonics
        [~, int] = integrate_fun_trimesh(elements.points', elements.otri,  squeeze(abs(squeeze(real(sum(H(k,:,:)))) - squeeze(real(sum(H(k-1,:,:))))).^2 )' );
        L2errors(i,k-1) = sqrt(int);    
    end
    U = squeeze(U(cN,:,:));
    numIterations(i) = cN;
    numDOF(i) = size(elements.tri,1);
    [~, int] = integrate_fun_trimesh(elements.points', elements.otri,  squeeze(abs(squeeze(real(sum(U(cN,:,:)))) - squeeze(real(sum(U(cN-1,:,:))))).^2 )' );
    error(i) = sqrt(int);
    [~, int] = integrate_fun_trimesh(elements.points', elements.otri,  squeeze(abs(squeeze(real(sum(U(cN,:,:))))).^2)' );
    L2Normt0(i) = sqrt(int);
    degen(i)= max(max(abs(f.* squeeze(abs(squeeze(real(sum(U(cN,:,:)))))))));
end