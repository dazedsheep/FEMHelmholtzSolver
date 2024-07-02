%% Simulation of the nonlinear periodic Westervelt equation with homogeneous nonlinearity parameter
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
sourceValueDomain = 4;

% point scatterers and their domain
values = [0];
refractionIndex = [1, 1];
% linear case
%values = [0, 0];
radii = [0.025];
centers = [0; 0];

diffusivity = 10^(-9);

minHarmonics = 12; % minimum number of harmonics
nHarmonics = 12; % maximum number of harmonics

gamma = 10^(-9);

beta = 1/speed_of_sound;

meshSize = 0.0005;

excitationPoints = [0;0.0];
% ultrasound pressure of the "point" source
pressure = 5*10^5;
excitationPointsSize = [0.001];


[elements] = initializeMultiLeveLSolver(meshSize, domain);

% plot the positions of excitation(s) and source(s)
objects = getGridPointsLE(elements, [centers excitationPoints], [radii excitationPointsSize]);
figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2), objects, 'facecolor', 'interp'); shading interp;
xlabel("x [m]");
ylabel("y [m]");

% construct non-linearity
f = constructF(elements, massDensity, speed_of_sound, refractionIndex, centers, radii, values, sourceValueDomain, true);

% construct all space dependent wave numbers for all harmonics
kappa = constructKappa(elements, diffusivity, speed_of_sound, omega, refractionIndex, centers, radii, values, nHarmonics);
%%
% build a gaussian source
% source = 1./(speed_of_sound.^2 + 1i .* omega .* diffusivity).*referencePressure.*gaussianSource(elements, excitationPoints, 0.6);
% build a point source (regularized dirac)

source = exp(1i.*omega.*pi/2).*pressure.*createPointSource(elements, excitationPoints, meshSize);  
excitation = zeros(size(elements.points,1),nHarmonics);
excitation(:,1) = source;

[cN, U, F] =  solveWesterveltMultiLevel(elements, omega, beta, gamma, kappa, excitation, f, nHarmonics, minHarmonics, false, 10^(-12));
H = U;
U = squeeze(U(cN,:,:));

%% Compute the solution
t_step = 1/(4/T);

%t_step = 0.0005;
% simulate 1000 steps
%tind = 0:t_step:10000*t_step;
tind = 0:0.01:0.02;
n_t = size(tind,2);

% compute only for n iterations
iter = 4;

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


%figure, plot(abs(p_L_2(2:(size(p_L_2,2))) - p_L_2(1:(size(p_L_2,2)-1))))