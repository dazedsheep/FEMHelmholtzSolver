%% Solver for the periodic Westervelt equation
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
sourceValueDomain = 5;

% point scatterers and their domain
values = [0, 0];
refractionIndex = [1, 1];
% linear case
%values = [0, 0];
radii = [0.1, 0.1];
centers = [0, 0; -0.3, 0.2];

diffusivity = 10^(-9);

minHarmonics = 12; % minimum number of harmonics
nHarmonics = 12; % maximum number of harmonics

% impdeance boundary conditions --> massDensity cancels
% higher frequencies are taken into account later
beta = 1/(speed_of_sound);
gamma = 10^(-9);

meshSize = 0.005;

% build a linear array
excitationPoints = [-2*lambda/8,-lambda/8,0,lambda/8,2*lambda/8;0.8,0.8,0.8,0.8,0.8];
% typical ultrasound pressure is 1MPa at a frequency of 1 MHz, lower
% frequency -> lower pressure!
pressure = 1*10^6;
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
%source = pressure.*gaussianSource(elements, excitationPoints, 0.6);
% build a point source (regularized dirac)

% we need to scale the reference pressure to the "point" source
source = pressure.*createPointSource(elements, excitationPoints, meshSize);  
excitation = zeros(size(elements.points,1),nHarmonics);
excitation(:,1) = source;

%excitation = 1i./(speed_of_sound.^2 + 1i .* (1:nHarmonics) .* omega .* diffusivity).*excitation;
[cN, U, F] = solveWesterveltMultiLevel(elements, omega, beta, gamma, kappa, excitation, f, nHarmonics, minHarmonics, false, 10^(-12));
H = U;
U = squeeze(U(cN,:,:));

%% linear case
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
sourceValueDomain = 0;

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
% typical ultrasound pressure is 1MPa at a frequency of 1 MHz, lower
% frequency -> lower pressure!
pressure = 1*10^6;
excitationPointsSize = [0.01,0.01,0.01,0.01,0.01];
excitationPower(1,1) = pressure;
excitationPower(1,2:nHarmonics) = 0;
nExcitationHarmonics = 1;

[elements] = initializeMultiLeveLSolver(meshSize, domain);

% plot the positions of excitation(s) and source(s)
objects = getGridPointsLE(elements, [centers excitationPoints], [radii excitationPointsSize]);
figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2), objects, 'facecolor', 'interp'); shading interp;
xlabel("x [m]");
ylabel("y [m]");

% construct nonlinearity
f = constructF(elements, massDensity, speed_of_sound, refractionIndex, centers, radii, values, sourceValueDomain, false);
% construct all space dependent wave numbers for all harmonics
kappa = constructKappa(elements, diffusivity, speed_of_sound, omega, refractionIndex, centers, radii, values, nHarmonics);

% build a gaussian source
%source = -speed_of_sound^2./(speed_of_sound.^2 + 1i .* omega .* diffusivity).*referencePressure.*gaussianSource(elements, excitationPoints, 0.6);
% build a point source (regularized dirac)

% we need to scale the reference pressure to the "point" source
source = omega^2./(speed_of_sound.^2 + 1i .* omega .* diffusivity).*createPointSource(elements, excitationPoints, meshSize).*pressure;  
excitation = zeros(size(elements.points,1),nHarmonics);
excitation(:,1) = source;

%excitation = 1i./(speed_of_sound.^2 + 1i .* (1:nHarmonics) .* omega .* diffusivity).*excitation;
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

%figure, plot(abs(p_L_2(2:(size(p_L_2,2))) - p_L_2(1:(size(p_L_2,2)-1))))

%% points along a line segment (cut through)
% lsegStart = [-1/2;-1/2];
% lsegEnd = [2/4;2/4];

% [v,peakidx] = max(P(1,1,:));
% 
% lsegStart = elements.points(peakidx,:)';
%excitationPoints(:,1) = [0;0.66];
lsegStart = [0;0.67];
lsegEnd = lsegStart + [0;-1.4];

npoints = 2000; 

unitVec = (lsegEnd - lsegStart)/norm((lsegEnd - lsegStart));
points = lsegStart + unitVec.*norm((lsegEnd - lsegStart),2)/npoints .* (0:(npoints-1));

idxList = zeros(npoints,1);

for j=1:size(points,2)
        [v,pcenterIdx] = min(sum((elements.points - points(:,j)').^2,2)); 
        idxList(j) = pcenterIdx;    
end

coords = [elements.points(idxList,1)';elements.points(idxList,2)'];

t = 1;

pPoint = squeeze(P(iter,t,idxList));

% find the points in space which are directly affected by nonlinearity
idxNonlinearity = find((sqrt(dot(points - centers(:,1) ,points - centers(:,1))) < radii(1))==1);

pointdist = sqrt(dot(points - lsegStart,points - lsegStart));
% smooth the data for plotting
smdata = smoothdata(pPoint, 'gaussian', 50);
figure, plot(pointdist, smdata)
hold on
plot(pointdist(idxNonlinearity), smdata(idxNonlinearity),'r')
xlabel('x');
ylabel('Pa');

%% points along a line segment (cut through)
% lsegStart = [-1/2;-1/2];
% lsegEnd = [2/4;2/4];

% lsegStart = excitationPoints(:,1);
% 
[v,peakidx] = max(P(1,1,:));
% 
lsegStart = elements.points(peakidx,:)';
%excitationPoints(:,1) = [0;0.66];
lsegStart = [0;0.67];
lsegEnd = lsegStart + [0;-1.4];

npoints = 4000; 

unitVec = (lsegEnd - lsegStart)/norm((lsegEnd - lsegStart));
points = lsegStart + unitVec.*norm((lsegEnd - lsegStart),2)/npoints .* (0:(npoints-1));

idxList = zeros(npoints,1);

for j=1:size(points,2)
        [v,pcenterIdx] = min(sum((elements.points - points(:,j)').^2,2)); 
        idxList(j) = pcenterIdx;    
end

coords = [elements.points(idxList,1)';elements.points(idxList,2)'];

t = 1;

pPoint = squeeze(P(iter,t,idxList));

% find the points in space which are directly affected by nonlinearity
idxNonlinearity = find((sqrt(dot(points - centers(:,1) ,points - centers(:,1))) < radii(1))==1);

pointdist = sqrt(dot(points - lsegStart,points - lsegStart));
% smooth the data for plotting
smdata = smoothdata(pPoint, 'gaussian', 50);
figure, plot(pointdist, smdata)
hold on
plot(pointdist(idxNonlinearity), smdata(idxNonlinearity),'r')
xlabel('x');
ylabel('Pa');
