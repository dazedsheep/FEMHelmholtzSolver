%%
% Impact of pollution effect on the convergence w.r.t. the element diameter
% kappa * h \leq 0.5 is a rule of thumb

clearvars

massDensity = 1000; %kg/m^3
speed_of_sound = 1450; % m/s

% signal period or center frequency
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
% prepare the source - take the smallest mesh and create our source
%[source_elements] = initializeMultiLeveLSolver(initialMeshsize/(meshdivfactor^(meshSteps-1)), domain);
% ultrasound pressure of the "point" source
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
    
    numDOF(i) = size(elements.triangles,2);

    [cN, U, F] =  solveWesterveltMultiLevel(elements, omega, beta, gamma, kappa, excitation, f, nHarmonics, minHarmonics, false, 0);
    numIterations(i) = cN;
    H = U;
    % calculate the difference of consecutive solution over the iterations
    for k=2:minHarmonics
        [~, int] = integrate_fun_trimesh(elements.points', elements.otri,  squeeze(abs(squeeze(real(sum(H(k,:,:)))) - squeeze(real(sum(H(k-1,:,:))))).^2 )' );
        L2errors(i,k-1) = sqrt(int);    
    end
    U = squeeze(U(cN,:,:));
    
  %  singleHelmholtz(i,:,:) = tri2grid(elements.opoints, elements.otri, real(squeeze(H(cN,1,:))), x, y); 
   
  %  sol(i,:,:) = tri2grid(elements.opoints, elements.otri, real(sum(U,1))', x, y);  
end

%%
% every point that is not in the domain is set to nan, set them 0 for the
% error calculations

sol(isnan(sol)) = 0;
singleHelmholtz(isnan(singleHelmholtz))= 0;

for i = 1:(meshSteps-1)
    error(i) = sqrt(sum(sum(squeeze(sol(i,:,:) - sol(i+1,:,:)).^2 .* (initialMeshsize/(meshdivfactor^(meshSteps+1))).^2)));
    nd(i) = numDOF(i+1)/numDOF(i);
    errorHelmholtz(i) = sqrt(sum(sum(squeeze(singleHelmholtz(i,:,:) - singleHelmholtz(i+1,:,:)).^2 .* (initialMeshsize/(meshdivfactor^(meshSteps+1))).^2)));
end

for i = 1: (length(error)-1)
    p(i) = log(error(i)/error(i+1))/log(2);
end

figure, plot(p);
%%

%%

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
%%
start = 1;

for i = start:(nHarmonics-1)
    ltwo(i - start + 1) = sqrt(sum(abs(num_sol(i+1,:) - num_sol(i,:)).^2));
    maxn(i - start + 1)  = max(abs(num_sol(i+1,:) - num_sol(i,:)));
    
end

figure, plot(log10(ltwo))
hold on
plot(log10(maxn))
legend("l2 norm", "max")
%% Next, we fix a certain error and increase the pressure and/or frequency
% pressure to show when the system blows up, frequency to show the impact on the grid 

clearvars

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

minHarmonics = 10; % minimum number of harmonics
nHarmonics = 40; % maximum number of harmonics

gamma = 1;

excitationPoints = [0;0];
excitationSize = 0.01;
% ultrasound pressure of the "point" source
pressure = 10;

% signal period or center frequency
T = 10^-3;

frequencySteps = 5;

initialMeshSize = 0.005;

% (power = 10, c = 50)
speed_of_sound = [50];

for j = 1:length(speed_of_sound)
    j
    for i = 1:frequencySteps
        
        beta = 1/speed_of_sound(j);
        
        omega = 2*pi* 1/T * i;
        
        meshSize(j,i) = initialMeshSize;
        
        while true      
            [elements] = initializeMultiLeveLSolver(meshSize(j,i), domain);
            %
            f = constructF(elements, massDensity, speed_of_sound(j), refractionIndex, centers, radii, values, sourceValueDomain, true);
    
            % construct all space dependent wave numbers for all harmonics
            kappa = constructKappa(elements, diffusivity, speed_of_sound(j), omega, refractionIndex, centers, radii, values, nHarmonics);
    
            % normalise the regularised dirac
            pointSource = createPointSource(elements, excitationPoints, excitationSize);
            [area, s] = integrate_fun_trimesh(elements.opoints, elements.otri, pointSource.');
            source = -exp(1i*pi/2).*pressure.*pointSource./s;  
            excitation = zeros(size(elements.points,1),nHarmonics);
            excitation(:,1) = source;
            numDOF(j,i) = size(elements.tri,1);
            [cN, U, F] =  solveWesterveltMultiLevel(elements, omega, beta, gamma, kappa, excitation, f, nHarmonics, minHarmonics, false, fixedError);
            numIterations(j,i) = cN;
            H = U;
            U = squeeze(U(cN,:,:));
            % compute the L^2 error
            [~, int] = integrate_fun_trimesh(elements.points', elements.otri,  squeeze(abs(squeeze(real(sum(H(cN,:,:)))) - squeeze(real(sum(H(cN-1,:,:))))).^2 )' );
            % we could also check whether cN < maximum number of iterations
            if sqrt(int) < fixedError
                break;
            end
            % if we did not reach the error level increase the resolution
            meshSize(j,i) = meshSize(i) / 2;
        end
        i
    end
    pressure = pressure * 2;
  
end

%% same for higher speed of sound

clearvars

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

minHarmonics = 10; % minimum number of harmonics
nHarmonics = 60; % maximum number of harmonics

gamma = 1;

excitationPoints = [0;0];
excitationSize = 0.01;
% ultrasound pressure of the "point" source
pressure = 10^6;

% signal period or center frequency
T = 10^-3;

frequencySteps = 5;

initialMeshSize = 0.005;

% (power = 10, c = 50)
speed_of_sound = [1450];

for j = 1:length(speed_of_sound)
    j
    for i = 1:frequencySteps
        
        beta = 1/speed_of_sound(j);
        
        omega = 2*pi* 1/T * i;
        
        meshSize(j,i) = initialMeshSize;
        
        while true      
            [elements] = initializeMultiLeveLSolver(meshSize(j,i), domain);
            %
            f = constructF(elements, massDensity, speed_of_sound(j), refractionIndex, centers, radii, values, sourceValueDomain, true);
    
            % construct all space dependent wave numbers for all harmonics
            kappa = constructKappa(elements, diffusivity, speed_of_sound(j), omega, refractionIndex, centers, radii, values, nHarmonics);
    
            % normalise the regularised dirac
            pointSource = createPointSource(elements, excitationPoints, excitationSize);
            [area, s] = integrate_fun_trimesh(elements.opoints, elements.otri, pointSource.');
            source = -exp(1i*pi/2).*pressure.*pointSource./s;  
            excitation = zeros(size(elements.points,1),nHarmonics);
            excitation(:,1) = source;
            numDOF(j,i) = size(elements.tri,1);
            [cN, U, F] =  solveWesterveltMultiLevel(elements, omega, beta, gamma, kappa, excitation, f, nHarmonics, minHarmonics, false, fixedError);
            H = U;
            U = squeeze(U(cN,:,:));
            numIterations(j,i) = cN;
            % compute the L^2 error
            [~, int] = integrate_fun_trimesh(elements.points', elements.otri,  squeeze(abs(squeeze(real(sum(H(cN,:,:)))) - squeeze(real(sum(H(cN-1,:,:))))).^2 )' );
            % we could also check whether cN < maximum number of iterations
            if sqrt(int) < fixedError
                break;
            end
            % if we did not reach the error level increase the resolution
            meshSize(j,i) = meshSize(i) / 2;
        end
        i
    end
    pressure = pressure * 2;
  
end
%%

freqvector = 1/T*(1:frequencySteps);
figure, plot(meshSize(1,:), freqvector);

%% compute the "solutions" and do a self convergence test (it suffices to do this for t = 0)

for i = 1:cN
    num_sol(i,:) = real(sum(H(i,:,:)));    
end

start = 1;

for i = start:(cN-1)
    [area, int] = integrate_fun_trimesh(elements.points', elements.otri,  squeeze((abs(num_sol(i+1,:) - num_sol(i,:)).^2)));
    ltwo(i - start + 1) = sqrt( int ) ;
    maxn(i - start + 1)  = max(abs(num_sol(i+1,:) - num_sol(i,:)));
end

figure, plot(log10(ltwo))
hold on
plot(log10(maxn))
legend("l2 norm", "max")

%% 
% run for certain amount of iterations and check self convergence (fixed
% mesh)

clearvars

fixedError = 10^-6; % L2 error of two consecutive iterations, we stop if we reach this

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

[cN, U, F] =  solveWesterveltMultiLevel(elements, omega, beta, gamma, kappa, excitation, f, nHarmonics, minHarmonics, false, 0);
H = U;
U = squeeze(U(cN,:,:));

%% compute the "solutions" and do a self convergence test (it suffices to do this for t = 0)

for i = 1:nHarmonics
    num_sol(i,:) = real(sum(H(i,:,:)));    
end

start = 1;

for i = start:(nHarmonics-1)
    [area, int] = integrate_fun_trimesh(elements.points', elements.otri,  squeeze((abs(num_sol(i+1,:) - num_sol(i,:)).^2)));
    ltwo(i - start + 1) = sqrt( int ) ;
    maxn(i - start + 1)  = max(abs(num_sol(i+1,:) - num_sol(i,:)));
    
end

figure, plot(log10(ltwo))
hold on
plot(log10(maxn))
legend("l2 norm", "max")

%%

clearvars

massDensity = 1000; %kg/m^3
speed_of_sound = 10; % m/s

% signal period or center frequency
T = 10^-2;
omega = 2*pi*1/T;
diffusivity = 10^(-9);

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

%%



sol(isnan(sol)) = 0;


for i = 1:(meshSteps-1)
    error(i) = sqrt(sum(sum(squeeze(sol(i,:,:) - sol(i+1,:,:)).^2.*(initialMeshsize/(2^meshSteps)))));


end

for i = 1: (length(error)-1)
    p(i) = log(error(i)/error(i+1))/log(2);
end

figure, plot(p);% the order should be approx 1 (L^2) for finer meshes
