%% This file simulates the study of the paper
% after computing the solutions for different inclusions and B/A values
% the figures are plotted
clearvars

massDensity = 1000; %kg/m^3
speed_of_sound = 344; % m/s (reference in air)
speed_of_sound_water = 1480; % m/s

% signal period or center frequency
T = 10^-4;
omega = 2*pi*1/T; % angular velocity

% wavelength
lambda = speed_of_sound/(1/T);

linArrayY = 0.0;
excitationPoints = [-4*lambda/4, -3*lambda/4, -2*lambda/4, -lambda/4, 0, lambda/4, 2*lambda/4, 3*lambda/4, 4*lambda/4;linArrayY,linArrayY,linArrayY,linArrayY,linArrayY,linArrayY,linArrayY,linArrayY,linArrayY];

% our domain
bcenter = [0,0];
brad = 1;
domain = [bcenter, brad];
% non linearity parameter of our domain (air = 1)
sourceValueDomain = 1;

% point scatterers and their domain
values = [5];
refractionIndex = [speed_of_sound/speed_of_sound_water];
radii = [0.35];
centers = [0; 0.37];

diffusivity = 10^(-9);

minHarmonics = 6; % minimum number of harmonics
nHarmonics = 6; % maximum number of harmonics

% impdeance boundary conditions --> massDensity cancels
% higher frequencies are taken into account later
beta = 1/(speed_of_sound);
gamma = 10^(-9);

meshSize = 0.002;
%linArrayY = 0.0;
% build a linear array
% excitationPoints = [-2*lambda/8,-lambda/8,0,lambda/8,2*lambda/8;linArrayY,linArrayY,linArrayY,linArrayY,linArrayY];
pressure = 4*10^4;
excitationPointsSize = [0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01];
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

% we need to scale the reference pressure to the "point" source
source = exp(1i.*omega.*pi/2).*pressure.*createPointSource(elements, excitationPoints, meshSize);  
excitation = zeros(size(elements.points,1),nHarmonics);
excitation(:,1) = source;

[cN, U, F] = solveWesterveltMultiLevel(elements, omega, beta, gamma, kappa, excitation, f, nHarmonics, minHarmonics, false, 10^(-12));
H = U;
U_no_phantoms = squeeze(U(cN,:,:));

P_no_phantoms = calcPressureProfile(omega, T, H, U_no_phantoms, cN);

clear H


% Phantom 1
% point scatterers and their domain
values = [5, 9];
refractionIndex = repmat(speed_of_sound/speed_of_sound_water, 1, 2);
radii = [0.35, 0.05];
centers = [0, -0.2; 0.37, 0.3];

% construct nonlinearity
f = constructF(elements, massDensity, speed_of_sound, refractionIndex, centers, radii, values, sourceValueDomain, true);
% construct all space dependent wave numbers for all harmonics
kappa = constructKappa(elements, diffusivity, speed_of_sound, omega, refractionIndex, centers, radii, values, nHarmonics);

% we need to scale the reference pressure to the "point" source
source = exp(1i.*omega.*pi/2).*pressure.*createPointSource(elements, excitationPoints, meshSize);  
excitation = zeros(size(elements.points,1),nHarmonics);
excitation(:,1) = source;

[cN, U, F] = solveWesterveltMultiLevel(elements, omega, beta, gamma, kappa, excitation, f, nHarmonics, minHarmonics, false, 10^(-12));
H = U;
U_one_phantom = squeeze(U(cN,:,:));

P_one_phantom = calcPressureProfile(omega, T, H, U_one_phantom, cN);

clear H

% Phantom 2

% point scatterers and their domain
values = [5, 10];
refractionIndex = repmat(speed_of_sound/speed_of_sound_water, 1, 2);

radii = [0.35, 0.04];
centers = [0, 0.2; 0.37, 0.3];

% construct nonlinearity
f = constructF(elements, massDensity, speed_of_sound, refractionIndex, centers, radii, values, sourceValueDomain, true);
% construct all space dependent wave numbers for all harmonics
kappa = constructKappa(elements, diffusivity, speed_of_sound, omega, refractionIndex, centers, radii, values, nHarmonics);

% we need to scale the reference pressure to the "point" source
source = exp(1i.*omega.*pi/2).*pressure.*createPointSource(elements, excitationPoints, meshSize);  
excitation = zeros(size(elements.points,1),nHarmonics);
excitation(:,1) = source;

[cN, U, F] = solveWesterveltMultiLevel(elements, omega, beta, gamma, kappa, excitation, f, nHarmonics, minHarmonics, false, 10^(-12));
H = U;
U_phantom_two = squeeze(U(cN,:,:));

P_phantom_two = calcPressureProfile(omega, T, H, U_phantom_two, cN);

clear H

% Phantom 3

% point scatterers and their domain
values = [5, 12];
refractionIndex = repmat(speed_of_sound/speed_of_sound_water, 1, 2);
radii = [0.35, 0.06];
centers = [0, 0; 0.37, 0.15];

% construct nonlinearity
f = constructF(elements, massDensity, speed_of_sound, refractionIndex, centers, radii, values, sourceValueDomain, true);
% construct all space dependent wave numbers for all harmonics
kappa = constructKappa(elements, diffusivity, speed_of_sound, omega, refractionIndex, centers, radii, values, nHarmonics);

% we need to scale the reference pressure to the "point" source
source = exp(1i.*omega.*pi/2).*pressure.*createPointSource(elements, excitationPoints, meshSize);  
excitation = zeros(size(elements.points,1),nHarmonics);
excitation(:,1) = source;

[cN, U, F] = solveWesterveltMultiLevel(elements, omega, beta, gamma, kappa, excitation, f, nHarmonics, minHarmonics, false, 10^(-12));
H = U;
U_phantom_three = squeeze(U(cN,:,:));

P_phantom_three = calcPressureProfile(omega, T, H, U_phantom_three, cN);

clear H


% 2 phantoms

% point scatterers and their domain
values = [5, 9, 10];
refractionIndex = repmat(speed_of_sound/speed_of_sound_water, 1, 3);

radii = [0.35, 0.05, 0.04];
centers = [0, -0.2, 0.2; 0.37, 0.3,0.3];

% construct nonlinearity
f = constructF(elements, massDensity, speed_of_sound, refractionIndex, centers, radii, values, sourceValueDomain, true);
% construct all space dependent wave numbers for all harmonics
kappa = constructKappa(elements, diffusivity, speed_of_sound, omega, refractionIndex, centers, radii, values, nHarmonics);


% we need to scale the reference pressure to the "point" source
source = exp(1i.*omega.*pi/2).*pressure.*createPointSource(elements, excitationPoints, meshSize);  
excitation = zeros(size(elements.points,1),nHarmonics);
excitation(:,1) = source;

[cN, U, F] = solveWesterveltMultiLevel(elements, omega, beta, gamma, kappa, excitation, f, nHarmonics, minHarmonics, false, 10^(-12));
H = U;
U_two_phantoms = squeeze(U(cN,:,:));

P_two_phantoms = calcPressureProfile(omega, T, H, U_two_phantoms, cN);

clear H

% 3 phantoms
% point scatterers and their domain
values = [5, 9, 10, 12];
refractionIndex = repmat(speed_of_sound/speed_of_sound_water, 1, 4);
radii = [0.35, 0.05, 0.04,0.06];
centers = [0, -0.2, 0.2, 0; 0.37, 0.3,0.3,0.15];

% construct nonlinearity
f = constructF(elements, massDensity, speed_of_sound, refractionIndex, centers, radii, values, sourceValueDomain, true);
% construct all space dependent wave numbers for all harmonics
kappa = constructKappa(elements, diffusivity, speed_of_sound, omega, refractionIndex, centers, radii, values, nHarmonics);

% we need to scale the reference pressure to the "point" source
source = exp(1i.*omega.*pi/2).*pressure.*createPointSource(elements, excitationPoints, meshSize);  
excitation = zeros(size(elements.points,1),nHarmonics);
excitation(:,1) = source;

[cN, U, F] = solveWesterveltMultiLevel(elements, omega, beta, gamma, kappa, excitation, f, nHarmonics, minHarmonics, false, 10^(-12));
H = U;
U_three_phantoms = squeeze(U(cN,:,:));

P_three_phantoms = calcPressureProfile(omega, T, H, U_three_phantoms, cN);

clear U H

%% Plot the boundary of our domain

boundaryValues = squeeze(P(1,1,elements.bedges(:,1)));
boundaryValuesFirstHarmonic = squeeze(lP(1,elements.bedges(:,1)));
figure, plot(abs(boundaryValues))
hold on
plot(abs(boundaryValuesFirstHarmonic))
hold off

figure, plot(boundaryValuesFirstHarmonic.' - boundaryValues);

scatter3(elements.points(elements.bedges(:,1),1), elements.points(elements.bedges(:,1),2), boundaryValues, 'filled');
hold on
scatter3(elements.points(elements.bedges(:,1),1), elements.points(elements.bedges(:,1),2), boundaryValuesFirstHarmonic, 'filled');

%% calc the point of the boundary of the object in our domain (water filled object)

% which of the objects is the domain of interest?
j = 1;

% get the elements for the object
objectVertices = [];
n = 1;
for i=1:size(elements.points,1)
    if norm(elements.points(i,:) - centers(:,j)',2) <= radii(j)
        objectVertices(n,:) = elements.points(i,:);
        objectVerticesIds(n,1) = i;
        n = n + 1;
    end
end

% calculate the objects' boundary
objectBoundary = boundary(objectVertices);
boundaryIndices = objectVerticesIds(objectBoundary,1);

% get the boundary values for the different solutions
boundaryValuesNoPhantoms = squeeze(P_no_phantoms(1,1,boundaryIndices));
boundaryValuesOnePhantom  = squeeze(P_one_phantom(1,1,boundaryIndices));
boundaryValuesTwoPhantoms  = squeeze(P_two_phantoms(1,1,boundaryIndices));
boundaryValuesThreePhantoms  = squeeze(P_three_phantoms(1,1,boundaryIndices));
boundaryValuesPhantomTwo = squeeze(P_phantom_two(1,1,boundaryIndices));
boundaryValuesPhantomThree = squeeze(P_phantom_three(1,1,boundaryIndices));

% find the boundary point nearest to the source
sourceLocation = [0,0];
dists = repmat(sourceLocation,size(elements.points(boundaryIndices,:),1),1)  - elements.points(boundaryIndices,:);
[dist, idx] = min(diag(dists*dists'));

% compute angles w.r.t. to source being at 0 rads
relativeBoundaryPoints = elements.points(boundaryIndices,:) - repmat(centers(:,1)',size(elements.points(boundaryIndices,:),1),1);

% get the angle and shift appropriately
angles = angle((relativeBoundaryPoints(:,1)+1i*relativeBoundaryPoints(:,2))*exp(1i*pi/2));

%%
figure, plot(abs(boundaryValuesNoPhantoms - boundaryValuesOnePhantom))
hold on
plot(abs(boundaryValuesNoPhantoms - boundaryValuesTwoPhantoms))
hold on
plot(abs(boundaryValuesNoPhantoms - boundaryValuesThreePhantoms))
hold on
xline(idx)

figure, scatter3(elements.points(boundaryIndices,1), elements.points(boundaryIndices,2), boundaryValuesNoPhantoms, 'filled');
xlabel("x")
ylabel("y")

% scatter differnce plots
figure, scatter3(elements.points(boundaryIndices,1), elements.points(boundaryIndices,2), boundaryValuesNoPhantoms - boundaryValuesOnePhantom, 'filled');
xlabel("x")
ylabel("y")

%% (difference) boundary values at t=0
indices = 1:(size(angles,1)-3); % this is due to the triangulation
b1 = boundaryValuesNoPhantoms - boundaryValuesOnePhantom;
b2 = boundaryValuesNoPhantoms - boundaryValuesPhantomTwo;
b3 = boundaryValuesNoPhantoms - boundaryValuesPhantomThree;
b4 = boundaryValuesNoPhantoms - boundaryValuesTwoPhantoms;
b5 = boundaryValuesNoPhantoms - boundaryValuesThreePhantoms;
figure
subplot(3,2,1)
plot(angles(indices), b1(indices))
xline(0)
xlabel('rad from source')
legend('$p_0(0,x) - p_1(0,x)$', 'Interpreter', 'latex')
xlim([-pi,pi])
ylabel('Pa')
subplot(3,2,3)
plot(angles(indices), b2(indices))
xline(0)
xlabel('rad from source')
ylabel('Pa')
legend('$p_0(0,x) - p_2(0,x)$', 'Interpreter', 'latex')
xlim([-pi,pi])

subplot(3,2,5)
plot(angles(indices), b3(indices))
xline(0)
xlabel('rad from source')
ylabel('Pa')
legend('$p_0(0,x) - p_3(0,x)$', 'Interpreter', 'latex')
xlim([-pi,pi])

subplot(3,2,2)
plot(angles(indices), boundaryValuesNoPhantoms(indices))
xline(0)
xlabel('rad from source')
ylabel('Pa')
legend('$p_0(0,x)$', 'Interpreter', 'latex')
xlim([-pi,pi])

subplot(3,2,4)
plot(angles(indices), b4(indices))
xline(0)
xlabel('rad from source')
ylabel('Pa')
legend('$p_0(0,x) - p_{1,2}(0,x)$', 'Interpreter', 'latex')
xlim([-pi,pi])

subplot(3,2,6)
plot(angles(indices), b5(indices))
xline(0)
xlabel('rad from source')
ylabel('Pa')
legend('$p_0(0,x) - p_{1,2,3}(0,x)$', 'Interpreter', 'latex')
xlim([-pi,pi])

%% auxilar plots
% select a point on the boundary and take look at the
% time behaviour

% sampling frequency
Fs = 10^7;
pointIdx = boundaryIndices(1157);
N = 4000;
pC_no_phantoms = zeros(1,N);
pC_two_phantoms = zeros(1,N); 
pC_three_phantoms = zeros(1,N);
for m=1:nHarmonics
    pC_no_phantoms = pC_no_phantoms + U_no_phantoms(m,pointIdx) .* exp(1i.*m.*omega.*(0:(N-1))*1/Fs);
    pC_three_phantoms = pC_three_phantoms + U_three_phantoms(m,pointIdx) .* exp(1i.*m.*omega.*(0:(N-1))*1/Fs);
    pC_two_phantoms = pC_two_phantoms + U_two_phantoms(m,pointIdx) .* exp(1i.*m.*omega.*(0:(N-1))*1/Fs);
end

time = (0:(N-1))*1/Fs;
timescale = 10^3;
figure
subplot(3,1,1)
plot(time*timescale,real(pC_no_phantoms))
ylabel("Acoustinc Pressure [Pa]");
subplot(3,1,2)
plot(time*timescale,real(pC_no_phantoms - pC_two_phantoms))
ylabel("Acoustinc Pressure [Pa]");
subplot(3,1,3)
plot(time*timescale,real(pC_no_phantoms - pC_three_phantoms))
ylabel("Acoustinc Pressure [Pa]");
xlabel("Time [ms]");


%% plot pressure profile at t=0
figure
subplot(1,2,1)
trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2), P_three_phantoms(1,1,:), 'facecolor', 'interp'); shading interp;
xlabel('x [m]')
ylabel('y [m]')
zlabel('Pa')
colormap('gray')
view(2)
subplot(1,2,2)
trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2), P_three_phantoms(1,1,:), 'facecolor', 'interp'); shading interp;
xlabel('x [m]')
ylabel('y [m]')
zlabel('Pa')
view([-37.5 30])
colormap('gray')
