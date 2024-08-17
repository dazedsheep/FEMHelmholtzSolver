%% This file simulates the first scenario of the paper including all plots
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
% nonlinearity parameter of our domain (water = 5)
sourceValueDomain = 5;

% point scatterers and their domain
values = [0];
refractionIndex = [1, 1];
% linear case
%values = [0, 0];
radii = [0.25];
centers = [0; 0];

diffusivity = 10^(-9);

minHarmonics = 5; % minimum number of harmonics
nHarmonics = 5; % maximum number of harmonics

gamma = 10^(-9);

beta = 1/speed_of_sound;

meshSize = 0.0005;

excitationPoints = [0.0,0.0;-0.2,0.2];
% ultrasound pressure of the "point" source
pressure = 3*10^7;
excitationPointsSize = [0.001;0.001];

[elements] = initializeMultiLeveLSolver(meshSize, domain);
%%
% construct non-linearity
f = constructF(elements, massDensity, speed_of_sound, refractionIndex, centers, radii, values, sourceValueDomain, true);

% construct all space dependent wave numbers for all harmonics
kappa = constructKappa(elements, diffusivity, speed_of_sound, omega, refractionIndex, centers, radii, values, nHarmonics);

source = exp(1i.*omega.*pi/2).*pressure.*createPointSourceOnBoundary(elements, excitationPoints, excitationPointsSize, meshSize);  
excitation = zeros(size(elements.points,1),nHarmonics);
excitation(:,1) = source;

[cN, U, F] =  solveWesterveltMultiLevelBoundaryExcitation(elements, omega, beta, gamma, kappa, excitation, f, nHarmonics, minHarmonics, false, 10^(-12));
H = U;
U = squeeze(U(cN,:,:));

P_excitation_on_the_boundary = calcPressureProfile(omega, T, H, U, cN);

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

%% Plot the frequency components 
point = [-0.12;0.0];

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