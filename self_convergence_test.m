%%
% Self-convergence w.r.t element diameter at fixed speed of sound
% and pressure, we should see p \approx 2, our method is of second order

clearvars

massDensity = 1000; %kg/m^3
speed_of_sound = 10; % m/s

% signal period or center frequency
T = 10^-2;
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

minHarmonics = 8; % minimum number of harmonics
nHarmonics = 8; % maximum number of harmonics

gamma = 10^(-9);

beta = 1/speed_of_sound;

meshSteps = 7;
initialMeshsize = 0.01;
x = (-0.1):(initialMeshsize/(2^meshSteps)):(0.1);
y = x;
for i = 1:meshSteps

    meshSize = initialMeshsize/2^(i-1);

    excitationPoints = [0;0];
    % ultrasound pressure of the "point" source
    pressure = 1*10^3;
    excitationPointsSize = [0.001];

    [elements] = initializeMultiLeveLSolver(meshSize, domain);


    f = constructF(elements, massDensity, speed_of_sound, refractionIndex, centers, radii, values, sourceValueDomain, true);

    % construct all space dependent wave numbers for all harmonics
    kappa = constructKappa(elements, diffusivity, speed_of_sound, omega, refractionIndex, centers, radii, values, nHarmonics);

    source = -exp(1i*pi/2).*pressure.*createPointSource(elements, excitationPoints, meshSize);  
    excitation = zeros(size(elements.points,1),nHarmonics);
    excitation(:,1) = source;

    [cN, U, F] =  solveWesterveltMultiLevel(elements, omega, beta, gamma, kappa, excitation, f, nHarmonics, minHarmonics, false, 10^(-12));
    H = U;
    U = squeeze(U(cN,:,:));

    sol(i,:,:) = tri2grid(elements.opoints, elements.otri, real(sum(U,1))', x, y);  
end

%
sol(isnan(sol)) = 0;
for i = 1:(meshSteps-1)
    error(i) = sqrt(sum(sum(squeeze(sol(i,:,:) - sol(i+1,:,:)).^2)));
end

for i = 1: (length(error)-1)
    p(i) = log(error(i)/error(i+1))/log(2);
end

figure, plot(p);% the order should be approx 2 (L^2) for finer meshes (as our method is of second order)

%% Next, we fix a certain error and increase the pressure and/or frequency
% pressure to show when the system blows up, frequency to show the impact on the grid 

clearvars

massDensity = 1000; %kg/m^3
speed_of_sound = 10; % m/s

% signal period or center frequency
T = 10^-2;
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

minHarmonics = 35; % minimum number of harmonics
nHarmonics = 35; % maximum number of harmonics

gamma = 10^(-9);

beta = 1/speed_of_sound;

meshSize = 0.0005;

excitationPoints = [0;0];
% ultrasound pressure of the "point" source
pressure = 1*10^3;
excitationPointsSize = [0.001];

[elements] = initializeMultiLeveLSolver(meshSize, domain);

%%
f = constructF(elements, massDensity, speed_of_sound, refractionIndex, centers, radii, values, sourceValueDomain, true);

% construct all space dependent wave numbers for all harmonics
kappa = constructKappa(elements, diffusivity, speed_of_sound, omega, refractionIndex, centers, radii, values, nHarmonics);

source = -exp(1i*pi/2).*pressure.*createPointSource(elements, excitationPoints, meshSize);  
excitation = zeros(size(elements.points,1),nHarmonics);
excitation(:,1) = source;

[cN, U, F] =  solveWesterveltMultiLevel(elements, omega, beta, gamma, kappa, excitation, f, nHarmonics, minHarmonics, false, 10^(-12));
H = U;
U = squeeze(U(cN,:,:));
%% compute the "solutions" and do a self convergence test (it suffices to do this for t = 0)

for i = 1:nHarmonics
    num_sol(i,:) = real(sum(H(i,:,:)));    
end

start = 1;

for i = start:(nHarmonics-1)
    ltwo(i - start + 1) = sqrt(sum(abs(num_sol(i+1,:) - num_sol(i,:)).^2));
    maxn(i - start + 1)  = max(abs(num_sol(i+1,:) - num_sol(i,:)));
    
end

figure, plot(log10(ltwo))
hold on
plot(log10(maxn))
legend("l2 norm", "max")