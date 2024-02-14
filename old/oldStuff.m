%% disc sources with multilevel scheme - with multiple sources near the boundary

clearvars

% excitation_points = [1/4,1/4;1/4,3/4];
% excitation_power_multiplier = 100;
% excitation_points_size = [0.03,0.03];

%excitation_points = [1/4;1/4];
%excitation_power_multiplier = 40;
%excitation_points_size = [0.07];

speed_of_sound = 1540; % m/s
%speed_of_sound = 15; % m/s

% signal period or center frequency
T = 10^-4;
w = 2*pi*1/T;

% our domain
bcenter = [1/2,1/2];
brad = 1/2;
domain = [bcenter, brad];

% place the excitation points uniformly near the boundary
num_points = 10;
bradmult = 0.9;

x = bradmult.*brad.*cos( (0:(num_points-1)) .* 2*pi / num_points) + bcenter(1);
y = bradmult.*brad.*sin( (0:(num_points-1)) .* 2*pi / num_points) + bcenter(1);
excitation_points = [x;y];
excitation_power_multiplier = -30;
excitation_points_size = ones(1,num_points).* 0.05;

% point scatterers and their domain
values = [4, 3];
radii = [0.15, 0.15];
centers = [3/4, 3/4; 1/2, 3/4];

% actually kappa = w/(sqrt(c^2 + i*w*b), b accounts for the diffusitivity
% of sound
b = 10;
kappa = w/sqrt(speed_of_sound^2 + 1i*w*b);
nHarmonics = 50;
gamma = 1;
beta = 1/speed_of_sound;
[boundaryIndices, elements, U, F, coupling] = solveForwardMultiLevel(speed_of_sound, w, gamma, beta, kappa, centers, radii, values, domain, excitation_points, excitation_points_size, excitation_power_multiplier, nHarmonics);
H = U;
U = squeeze(U(nHarmonics,:,:));

% difference between harmonics figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2), abs(real(H(6,3,:))-real(H(14,3,:))), 'facecolor', 'interp'); shading interp;

%% disc sources with multilevel scheme

clearvars

% excitation_points = [1/4,1/4;1/4,3/4];
% excitation_power_multiplier = 100;
% excitation_points_size = [0.03,0.03];

excitation_points = [-1/2;-1/2];
excitation_power_multiplier = -30;
excitation_points_size = [0.07];

speed_of_sound = 1540; % m/s
%speed_of_sound = 15; % m/s

% signal period or center frequency
T = 10^-4;
w = 2*pi*1/T;

% our domain
bcenter = [0,0];
brad = 1;
domain = [bcenter, brad];

% point scatterers and their domain
values = [4, 0];
radii = [0.2, 0.1];
centers = [0, 1/2; 0, -1/4];

% actually kappa = w/(sqrt(c^2 + i*w*b), b accounts for the diffusitivity
% of sound
b = 0.0001;
kappa = w/sqrt(speed_of_sound^2 + 1i*w*b);
minHarmonics = 10; % minimum number of harmonics
nHarmonics = 12; % maximum number of harmonics
gamma = 1;
beta = 1/speed_of_sound;
meshSize = 0.005;
[cN, boundaryIndices, elements, U, F, coupling] = solveForwardMultiLevelVectorized(w, gamma, beta, kappa, centers, radii, values, domain, excitation_points, excitation_points_size, excitation_power_multiplier, nHarmonics, minHarmonics, 10^(-12), meshSize);
H = U;
U = squeeze(U(cN,:,:));

% difference between harmonics figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2), abs(real(H(6,3,:))-real(H(14,3,:))), 'facecolor', 'interp'); shading interp;

%% disc sources with multilevel scheme - High frequency

clearvars

% excitation_points = [1/4,1/4;1/4,3/4];
% excitation_power_multiplier = 100;
% excitation_points_size = [0.03,0.03];

excitation_points = [-1/8;0];
excitation_power_multiplier = -40;
excitation_points_size = [0.07];

speed_of_sound = 1540; % m/s
%speed_of_sound = 15; % m/s

% signal period or center frequency
T = 10^-5;
w = 2*pi*1/T;

% our domain
bcenter = [0,0];
brad = 1/4;
domain = [bcenter, brad];

% point scatterers and their domain
values = [9, 0];

% linear case
%values = [0, 0];
radii = [0.02, 0.1];
centers = [0, 1/2; 0, -1/4];

% actually kappa = w/(sqrt(c^2 + i*w*b), b accounts for the diffusitivity
% of sound
b = 0.0001;
kappa = w/sqrt(speed_of_sound^2 + 1i*w*b);
minHarmonics = 4; % minimum number of harmonics
nHarmonics = 4; % maximum number of harmonics
gamma = 1;
beta = 1/speed_of_sound;
meshSize = 0.005;
[cN, boundaryIndices, elements, U, F, coupling] = solveForwardMultiLevelVectorized(w, gamma, beta, kappa, centers, radii, values, domain, excitation_points, excitation_points_size, excitation_power_multiplier, nHarmonics, minHarmonics, 10^(-12), meshSize);
H = U;
U = squeeze(U(cN,:,:));

% difference between harmonics figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2), abs(real(H(6,3,:))-real(H(14,3,:))), 'facecolor', 'interp'); shading interp;




%% disc sources with multilevel scheme - low frequency
clearvars
% excitation_points = [1/4,1/4;1/4,3/4];
% excitation_power_multiplier = 100;
% excitation_points_size = [0.03,0.03];

excitation_points = [-7;0];
excitation_power_multiplier = -40;
excitation_points_size = [0.02];

speed_of_sound = 1540; % m/s
%speed_of_sound = 15; % m/s

% signal period or center frequency
T = 10^-3;
w = 2*pi*1/T;

% our domain
bcenter = [0,0];
brad = 1;
domain = [bcenter, brad];

% point scatterers and their domain
%values = [9, 0];

% linear case
values = [4, 0];
radii = [0.1, 0.1];
centers = [0, 1/2; 0, -1/4];

% actually kappa = w/(sqrt(c^2 + i*w*b), b accounts for the diffusitivity
% of sound
b = 0.0001;
kappa = w/sqrt(speed_of_sound^2 + 1i*w*b);
minHarmonics = 20; % minimum number of harmonics
nHarmonics = 25; % maximum number of harmonics
gamma = 1;
beta = 1/speed_of_sound;
meshSize = 0.005;  % this could take long :)
[cN, boundaryIndices, elements, U, F, coupling] = solveForwardMultiLevelVectorized(w, gamma, beta, kappa, centers, radii, values, domain, excitation_points, excitation_points_size, excitation_power_multiplier, nHarmonics, minHarmonics, 10^(-12), meshSize);
H = U;
U = squeeze(U(cN,:,:));

% difference between harmonics figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2), abs(real(H(6,3,:))-real(H(14,3,:))), 'facecolor', 'interp'); shading interp;

%% point sources
clearvars
% speef of sounde
speed_of_sound = 1540; % m/s


% signal period or cener frequency
T = 10^-4;
w = 2*pi*1/T;

% our domain
bcenter = [1/2,1/2];
brad = 1/2;
domain = [bcenter, brad];


% point scatterers and their domain
values = [20000, 22000];
radii = [0, 0];
centers = [3/4, 3/4; 1/4, 3/4];

% actually kappa = w/(sqrt(c^2 + i*w*b), b accounts for the diffusitivity
% of sound
b = 1;
kappa = w/sqrt(speed_of_sound^2 + 1i*w*b);
gamma = 1;
nHarmonics = 7;
[boundaryIndices, elements, U, F, coupling] = solveForward(speed_of_sound, w, gamma, kappa, centers, radii, values, domain, nHarmonics);

%% disc sources
clearvars
speed_of_sound = 1540; % m/s
%speed_of_sound = 15; % m/s

% signal period or center frequency
T = 10^-4;
w = 2*pi*1/T;

% our domain
bcenter = [1/2,1/2];
brad = 1/2;
domain = [bcenter, brad];


% point scatterers and their domain
values = [300, 270];
radii = [0.05, 0.04];
centers = [3/4, 3/4; 1/4, 3/4];

% actually kappa = w/(sqrt(c^2 + i*w*b), b accounts for the diffusitivity
% of sound, we set b=0, neglecting the diffusitivity of sound
% (justification: \rho_0 is big and therefore b is very very very small)
kappa = w/speed_of_sound;

nHarmonics = 7;
gamma = 1;

[boundaryIndices, elements, U, F, coupling] = solveForward(speed_of_sound, w, gamma, kappa, centers, radii, values, domain, nHarmonics);

%% disc sources with fixed point scheme

clearvars
speed_of_sound = 1540; % m/s
%speed_of_sound = 15; % m/s

% signal period or center frequency
T = 10^-4;
w = 2*pi*1/T;

% our domain
bcenter = [1/2,1/2];
brad = 1/2;
domain = [bcenter, brad];


% point scatterers and their domain
values = [4, 2];
radii = [0.05, 0.04];
centers = [3/4, 3/4; 1/4, 3/4];

% actually kappa = w/(sqrt(c^2 + i*w*b), b accounts for the diffusitivity
% of sound, we set b=0, neglecting the diffusitivity of sound
% (justification: \rho_0 is big and therefore b is very very very small)
kappa = w/speed_of_sound;

iterations = 10;
gamma = 1;

[boundaryIndices, elements, U, F, coupling] = solveForwardFixedPoint(speed_of_sound, w, gamma, kappa, centers, radii, values, domain, iterations);
H = U;
U = squeeze(U(iterations,:,:));



%% plain gradient descent for identifying sources and their values
m = 2;
sourcePoints = [];
sourceValues = [];

f = zeros(size(elements.points,1),1);
y = U(2,elements.bedges(:,1));
g = f;
h = f;
f(elements.bedges(:,1)) = -1i.*m.*kappa.*y;
% K*(-y) 
%te = solveHelmholtzVectorizedTmp(elements, m*w, m*kappa, -1/speed_of_sound, zeros(size(elements.points,1),1), f, g, size(elements.points,1));

%figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2), real(te), 'facecolor', 'interp'); shading interp;
%figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2), real(U(2,:)), 'facecolor', 'interp'); shading interp;


for j = 1:2
    
    te = solveHelmholtzVectorizedTmp(elements, m*w, m*kappa, -1/speed_of_sound, zeros(size(elements.points,1),1), f, g, size(elements.points,1));

    [v, pc] = max(abs(te));
    center = elements.points(pc,:)';
    sourcePoints = [sourcePoints center];
    % try to maximize the centers source value ||Kx - y||^2, where x are our
    % source values
    v = 1;
    F = constructF(elements, sourcePoints, zeros(size(sourcePoints,2),1), [sourceValues v]);
    v = [sourceValues v];
    fun = @(x) solveHelmholtzVectorizedTmp(elements, m*w, m*kappa, 1/speed_of_sound, constructF(elements, sourcePoints, zeros(size(sourcePoints,2),1), [x]).*m^2.*kappa^2.*coupling(m,:), g, g, size(elements.points,1));
    sigma = 10^-4;
    for k = 1:200
        lRate = 1000./min(abs(coupling(m,:)));
        % gradient
        h = 0.001;
        x = v;
        for i=1:size(x,2)
            e_i = zeros(1,size(x,2));
            e_i(i) = 1;
            y0 = fun(x-h*e_i);
            z = (y0(elements.bedges(:,1)).' - y);
            val0 = z*z';
    
    %         y0 = fun(x);
    %         z = (y0(elements.bedges(:,1)).' - y);
    %         val0 = z*z';
            y1 = fun(x+h*e_i);
            z = (y1(elements.bedges(:,1)).' - y);
            val1 = z*z';
            df(i,k) = 1/(2*h) *(val1 - val0);  
        end
        
        direction = -df(:,k);
        y_current = fun(x);
        loss = (y_current(elements.bedges(:,1)).' - y)*(y_current(elements.bedges(:,1)).' - y)';
        
        % line search / amijo conditions
        while 1
            y_new = fun(x + lRate*direction');
            loss_new = (y_new(elements.bedges(:,1)).' - y)*(y_new(elements.bedges(:,1)).' - y)';
            if (loss_new < loss + sigma*lRate*(direction'*direction))
                break;
            end
            lRate = lRate/10;
        end
    
        v = x - lRate*df(:,k)';
    end
    
    sourceValues = v;
    f(elements.bedges(:,1)) = 1i.*m.*kappa.*(y_new(elements.bedges(:,1)).' - y);


end
% uu = solveHelmholtzVectorizedTmp(elements, m*w, m*kappa, 1/speed_of_sound, F, g, g, size(elements.points,1));
% f(elements.bedges(:,1)) = -1i.*m.*kappa.*(uu(elements.bedges(:,1)).' - y);
% te2 = solveHelmholtzVectorizedTmp(elements, m*w, m*kappa, -1/speed_of_sound, zeros(size(elements.points,1),1), f, g, size(elements.points,1));
