clear all

% center frequency of ultrasound
speed_of_sound = 1540; % m/s

% 20 kHz
% f = 20*10^3;
T = 10^-4;
w = 2*pi*1/T;

% our domain
bcenter = [1/2,1/2];
brad = 1/2;
domain = [bcenter, brad];


% point scatterers and their domain
values = [20, 30];
radii = [1/10, 1/20];
centers = [1/4, 3/4; 1/4, 3/4];

% values = [5];
% radii = [1/10];
% centers = [1/4; 1/4];

% actually kappa = w/(sqrt(c^2 + i*w*b), b accounts for the diffusitivity
% of sound, we set b=0, neglecting the diffusitivity of sound
kappa = w/speed_of_sound;

nHarmonics = 4;
[boundaryIndices, elements, U] = generateBoundaryData(speed_of_sound, w, kappa, centers, radii, values, domain, nHarmonics);
%%

t_step = 0.005;
n_t = floor(2*pi/t_step);

tind = 0:t_step:2*pi;

% construct p(t,x)
z = repmat(exp(-1i.*w.*T.*tind)', 1, size(U,2));
pC = zeros(size(z,1), size(z,2));

for k=1:nHarmonics
    pC = pC + repmat(U(k,:),size(z,1),1).*repmat(exp(-1i.*k.*w.*tind.*T)', 1, size(U,2));
end

p = real(pC);

%% let's take a look at p(t,x); playback the periodic solution
figure;
for i=1:n_t
    trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2), p(i,:), 'facecolor', 'interp'); shading interp;
    view([0 90]);
    drawnow
    pause(0.1)
end
