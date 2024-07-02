%% Solver for the periodic Westervelt equation
clearvars

massDensity = 1000; %kg/m^3
speed_of_sound = 344; % m/s (reference in water)

% signal period or center frequency
T = 10^-4;
omega = 2*pi*1/T; % angular velocity

% wavelength
lambda = speed_of_sound/(1/T);

linArrayY = 0.0;
excitationPoints = [-4*lambda/4, -3*lambda/4, -2*lambda/4,-lambda/4, 0, lambda/4, 2*lambda/4, 3*lambda/4, 4*lambda/4;linArrayY,linArrayY,linArrayY,linArrayY,linArrayY,linArrayY,linArrayY,linArrayY,linArrayY];
%excitationPoints = [0;0];

speed_of_sound = 344;


% our domain
bcenter = [0,0];
brad = 1;
domain = [bcenter, brad];
% non linearity parameter of our domain (air = 1)
sourceValueDomain = 1;

% point scatterers and their domain
values = [7, 9, 10, 15];
refractionIndex = [1/4.3023, 1/4.3023, 1/4.3023, 1/4.3023];
% linear case
%values = [0, 0];
radii = [0.35, 0.05, 0.04,0.06];
centers = [0, -0.2, 0.2, 0; 0.37, 0.3,0.3,0.15];

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

% build a gaussian source
%source = pressure.*gaussianSource(elements, excitationPoints, 0.6);
% build a point source (regularized dirac)

% we need to scale the reference pressure to the "point" source
source = exp(1i.*omega.*pi/2).*pressure.*createPointSource(elements, excitationPoints, meshSize);  
excitation = zeros(size(elements.points,1),nHarmonics);
excitation(:,1) = source;

[cN, U, F] = solveWesterveltMultiLevel(elements, omega, beta, gamma, kappa, excitation, f, nHarmonics, minHarmonics, false, 10^(-12));
H = U;
U = squeeze(U(cN,:,:));

%%
% calc the solution(s) for the multi level harmonic approximation
t_step = 1/(4/T);
%t_step = 0.0005;
% simulate 1000 steps
%tind = 0:t_step:10000*t_step;
tind = 0:0.01:0.02;
n_t = size(tind,2);

% compute only for n iterations
iter = 1;

% construct p(t,x)
z = repmat(exp(-1i.*omega.*T.*tind)', 1, size(U,2));
P = zeros(iter, size(z,1), size(z,2));
o = 0;
for j=(min(size(H,1),cN) - iter):min(size(H,1),cN)
    o = o + 1;
    rpC = zeros(size(z,1), size(z,2));
    for k=1:j
        rpC = rpC + squeeze(repmat(H(j,k,:),size(z,1),1)).*repmat(exp(1i.*k.*omega.*tind.*T).', 1, size(U,2));
    end
    p_L_2(o) = sqrt(sum(sum(real(rpC).^2, 2).^2,1));
    p_L_2_t(o) = sqrt(sum(real(rpC(1,:)).^2, 2));
    P(o,:,:) = real(rpC);
   
end

% in this case this is not the linear solution as we use convection
lP = real(squeeze(repmat(H(1,1,:),size(z,1),1)).*repmat(exp(1i.*omega.*tind.*T).', 1, size(U,2)));

%% pick a point and look at the frequency domain
% therefore we sample at 2*f_max which is in our case 2*nHarmonics*1/T
point = [0.0;0.1];

[v,idx] = min(sum((elements.points - point(:)').^2,2)); 

node = [elements.points(idx,1);elements.points(idx,2)];

% sampling frequency in time
Fs = 1/T * 2 * (nHarmonics+1);

N = 2000;
pC = zeros(1,N);
for m=1:nHarmonics
    pC = pC + U(m,idx) .* exp(1i.*m.*omega.*(0:(N-1))*1/Fs);
end

% time = (0:(N-1))*1/Fs;
% timescale = 10^3;
% figure, plot(time*timescale,real(pC))
% ylabel("Acoustinc Pressure [Pa]");
% xlabel("Time [ms]");

MaxBins = 5;
P0 = max(pressure,max(P(1,1,:)));


window = hanning(N);
freq =  Fs/N*(0:(N/2));
y = abs(fft(window'.*real(pC)))/N;
y1 = y(1:N/2+1);
y1(2:end-1) = 2*y1(2:end-1);
shiftedTF = fftshift(fft(real(pC)))/N;
TFdB = 10*log10(y1/P0);
fscaling = 10^3;
M = min(MaxBins, N/2 + 1);
xaxis = freq./fscaling;
figure, plot(xaxis, TFdB(1:N/2+1))
xlabel("Frequency [kHz]")
ylabel("P/P0 [dB]")
%%
Fs = 1/T * 2 * (nHarmonics+1)*5;

N = 2000;
pC = zeros(1,N);
for m=1:nHarmonics
    pC = pC + U(m,idx) .* exp(1i.*m.*omega.*(0:(N-1))*1/Fs);
end
lpC = H(1,1,idx).* exp(1i.*omega.*(0:(N-1))*1/Fs);

time = (0:(N-1))*1/Fs;
timescale = 10^3;
figure, plot(time*timescale,real(pC))
hold on
plot(time*timescale,real(lpC), "LineStyle","--")
ylabel("Acoustinc Pressure [Pa]");
xlabel("Time [ms]");

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

%% Plot the boundary of the object in our domain (water filled object)

%which of the objects is the domain of interest?
j = 1;

objectVertices = [];
n = 1;
for i=1:size(elements.points,1)
    if norm(elements.points(i,:) - centers(:,j)',2) < radii(j)
        objectVertices(n,:) = elements.points(i,:);
        objectVerticesIds(n,1) = i;
        n = n + 1;
    end
end

objectBoundary = boundary(objectVertices);
boundaryIndices = objectVerticesIds(objectBoundary,1);
boundaryValues = squeeze(P(1,1,boundaryIndices));
boundaryValuesFirstHarmonic = squeeze(lP(1,boundaryIndices));
figure, plot(abs(boundaryValues))
hold on
plot(abs(boundaryValuesFirstHarmonic))
hold off

figure, scatter3(elements.points(boundaryIndices,1), elements.points(boundaryIndices,2), boundaryValues, 'filled');
xlabel("x")
ylabel("y")
