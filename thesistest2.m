%% Solver for the periodic Westervelt equation
clearvars

massDensity = 1000; %kg/m^3
speed_of_sound = 1480; % m/s

% signal period or center frequency
T = 10^-5;
omega = 2*pi*1/T;

% our domain
bcenter = [0,0];
brad = 0.2;
domain = [bcenter, brad];
% non linearity parameter of our domain (water = 5)
sourceValueDomain = 0;

% point scatterers and their domain
values = [6, 10, 5, 13];
refractionIndex = [1,1,1,1];
radii = [0.02, 0.01, 0.015, 0.015];
centers = [-0.08, 0.08, 0, 0; 0, 0, -0.1, 0.05];

diffusivity = 10^(-9);

minHarmonics = 7; % minimum number of harmonics
nHarmonics = 7; % maximum number of harmonics

% impdeance boundary conditions --> massDensity cancels
% higher frequencies are taken into account later

gamma = 10^(-9);

beta = 1/speed_of_sound;

meshSize = 0.0005;

excitationPoints = [0;0];
% typical ultrasound pressure is 1MPa
pressure = 10^6;
excitationPointsSize = [0.001];

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
%source = 1./(speed_of_sound.^2 + 1i .* omega .* diffusivity).*referencePressure.*gaussianSource(elements, excitationPoints, 0.6);
% build a point source (regularized dirac)

% we need to scale the reference pressure to the "point" source/linear
% array
source = exp(1i.*pi).*pressure.*createPointSource(elements, excitationPoints, meshSize);
excitation = zeros(size(elements.points,1),nHarmonics);
excitation(:,1) = source;

[cN, U, F] =  solveWesterveltMultiLevel(elements, omega, beta, gamma, kappa, excitation, f, nHarmonics, minHarmonics, false, 10^(-12));
H = U;
U = squeeze(U(cN,:,:));


%% points along a line segment (cut through)
% lsegStart = [-1/2;-1/2];
% lsegEnd = [2/4;2/4];
% due to the mesh, the peak  of the source may be a little off
[v,peakidx] = max(P(1,1,:));

%lsegStart = elements.points(peakidx,:)';
lsegStart = [0;0];
lsegEnd = [0;0.2];

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
lpPoint = squeeze(lP(t, idxList));
% find the points in space which are directly affected by nonlinearity
idxNonlinearity = find((sqrt(dot(points - centers(:,1) ,points - centers(:,1))) < radii(1))==1);

pointdist = sqrt(dot(points - excitationPoints,points - excitationPoints));
% smooth the data for plotting
smdata = smoothdata(pPoint, 'gaussian', 50);
lpsmdata = smoothdata(lpPoint, 'gaussian', 50);
figure, plot(pointdist, smdata,"LineWidth",1)
hold on
plot(pointdist, lpsmdata, "LineStyle", "--", "Color", "black","LineWidth",1);
%plot(pointdist(idxNonlinearity), smdata(idxNonlinearity),'r')
xlabel('x [m]');
ylabel('Acoustic Pressure [Pa]');

%% Plot the boundary of our domain

boundaryValues = squeeze(P(1,1,elements.bedges(:,1)));
boundaryValuesLinear = squeeze(lP(1,elements.bedges(:,1)));
figure, plot(abs(boundaryValues))
hold on
plot(abs(boundaryValuesLinear))
hold off

figure, plot(boundaryValuesLinear.' - boundaryValues);

scatter3(elements.points(elements.bedges(:,1),1), elements.points(elements.bedges(:,1),2), boundaryValues, 'filled');
hold on
scatter3(elements.points(elements.bedges(:,1),1), elements.points(elements.bedges(:,1),2), boundaryValuesLinear, 'filled');
