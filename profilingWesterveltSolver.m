%% disc sources with multilevel scheme - with variable speed of sound

clearvars

massDensity = 1000; %kg/m^3

% excitation_points = [1/4,1/4;1/4,3/4];
% excitation_power_multiplier = 100;
% excitation_points_size = [0.03,0.03];

excitation_points = [-0.12;-0.12];
% typical ultrasound pressure is 1MPa
excitation_power_multiplier = -1*10^6;
excitation_points_size = [0.07];

speed_of_sound = 1540; % m/s

% signal period or center frequency
T = 1/5*10^-4;
w = 2*pi*1/T;

% our domain
bcenter = [0,0];
brad = 1/3;
domain = [bcenter, brad];

% point scatterers and their domain
values = [3.5, 0];
refractionIndex = [1.1, 1];
% linear case
%values = [0, 0];
radii = [0.05, 0.1];
centers = [0, 1/2; 0, -1/4];

% actually kappa = w/(sqrt(c^2 + i*w*b), b accounts for the diffusitivity
% of sound
b = 10^(-9);

kappa = w./sqrt(speed_of_sound.^2 + 1i*w*b);
minHarmonics = 6; % minimum number of harmonics
nHarmonics = 6; % maximum number of harmonics
gamma = 1;
beta = 1/speed_of_sound;
meshSize = 0.0005;
[cN, boundaryIndices, elements, U, F, coupling] = solveForwardMultiLevelCVectorized(w, gamma, beta, kappa, centers, radii, values, domain, excitation_points, excitation_points_size, excitation_power_multiplier, refractionIndex, speed_of_sound, b, nHarmonics, minHarmonics, 10^(-12),massDensity, meshSize);
H = U;
U = squeeze(U(cN,:,:));