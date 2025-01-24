%% This file simulates the intravascular study of the paper
% after computing the solutions for different inclusions and B/A values
% the figures are plotted
clearvars

massDensity = 1000; %kg/m^3
speed_of_sound = 344; % m/s (reference in air)
speed_of_sound_water = 1480; % m/s

% signal period or center frequency
T = 10^-6;
omega = 2*pi*1/T; % angular velocity

% wavelength
lambda = speed_of_sound/(1/T);

linArrayY = 0.0;
%excitationPoints = [-4*lambda/4, -3*lambda/4, -2*lambda/4, -lambda/4, 0, lambda/4, 2*lambda/4, 3*lambda/4, 4*lambda/4;linArrayY,linArrayY,linArrayY,linArrayY,linArrayY,linArrayY,linArrayY,linArrayY,linArrayY];
excitationPoints = [0;0];
% our domain
bcenter = [0,0];
brad = 0.01;
domain = [bcenter, brad];
% non linearity parameter of our domain (air = 1)
sourceValueDomain = 5;

% point scatterers and their domain
values = [5, 6, 7, 8, 11, 9];
refractionIndex = ones(1,6).* speed_of_sound_water/speed_of_sound_water;
radii = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001];
centers = [0 0 0.005 -0.005 0.005 -0.005; 0.005 -0.005 0 0 0.005 -0.005];

diffusivity = 10^(-9);

minHarmonics = 3; % minimum number of harmonics
nHarmonics = 3; % maximum number of harmonics

% impdeance boundary conditions --> massDensity cancels
% higher frequencies are taken into account later
beta = 1/(speed_of_sound_water);
gamma = 10^(-9);

meshSize = 0.00002;
%linArrayY = 0.0;
% build a linear array
% excitationPoints = [-2*lambda/8,-lambda/8,0,lambda/8,2*lambda/8;linArrayY,linArrayY,linArrayY,linArrayY,linArrayY];
pressure = 2*10^3;
excitationPointsSize = [0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001];
excitationPower(1,1) = pressure;
excitationPower(1,2:nHarmonics) = 0;

[elements] = initializeMultiLeveLSolver(meshSize, domain);
%%
% plot the positions of excitation(s) and source(s)
objects = getGridPointsLE(elements, [centers excitationPoints], [radii excitationPointsSize]);
figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2), objects, 'facecolor', 'interp'); shading interp;
xlabel("x [m]");
ylabel("y [m]");
%%
% construct nonlinearity
f = constructF(elements, massDensity, speed_of_sound, refractionIndex, centers, radii, values, sourceValueDomain, true);
% construct all space dependent wave numbers for all harmonics
kappa = constructKappa(elements, diffusivity, speed_of_sound, omega, refractionIndex, centers, radii, values, nHarmonics);

source = exp(1i.*omega.*-pi/2).*pressure.*createPointSource(elements, excitationPoints, meshSize);  
excitation = zeros(size(elements.points,1),nHarmonics);
excitation(:,1) = source;

[cN, U, F] =  solveWesterveltMultiLevel(elements, omega, beta, gamma, kappa, excitation, f, nHarmonics, minHarmonics, false, 10^(-12));
H = U;
U = squeeze(U(cN,:,:));

%% Compute the solution
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

% Calculate the solution to the linear equation for comparison
lP = real(squeeze(repmat(H(1,1,:),size(z,1),1)).*repmat(exp(1i.*omega.*tind.*T).', 1, size(U,2)));

%% Plot points along a line segment (cut through)

lsegStart = [-0.01;0];
lsegEnd = [0.01;0.0];
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

pPoint = squeeze(P(iter, t, idxList));
lpPoint = squeeze(lP(t, idxList));
% find the points in space which are directly affected by nonlinearity
idxNonlinearity = find((sqrt(dot(points - centers(:,1) ,points - centers(:,1))) < radii(1))==1);

pointdist = sqrt(dot(points - lsegStart,points - lsegStart));
% smooth the data for plotting
smdata = smoothdata(pPoint, 'gaussian', 50);
lpsmdata = smoothdata(lpPoint, 'gaussian', 50);
figure, plot(pointdist, smdata,"LineWidth",1)
hold on
plot(pointdist, lpsmdata, "LineStyle", "--", "Color", "black","LineWidth",1);
xlabel('x [m]');
ylabel('Acoustic Pressure [Pa]');

% show also the difference between linear and nonlinear solution
figure, plot(pointdist, smdata - lpsmdata);

%% Plot the frequency components 
point = [0.0;0.0];

[v,idx] = min(sum((elements.points - point(:)').^2,2)); 

node = [elements.points(idx,1);elements.points(idx,2)];

% sampling frequency in time
Fs = 1/T * 2 * (nHarmonics+1);

N = 2000;
pC = zeros(1,N);
for m=1:nHarmonics
    pC = pC + U(m,idx) .* exp(1i.*m.*omega.*(0:(N-1))*1/Fs);
end

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
