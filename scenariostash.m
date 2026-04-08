%% this works quite well, there is just not enough power in the second measurement (higher frequency)


% specify our reference states
f1 = 60;    % Hz
f2 = 123;    % Hz
f3 = f1;    % frequency of third reference state = frequency of first reference state
omega1 = 2*pi*f1;
omega2 = 2*pi*f2;
omega3 = 2*pi*f3;
u3Amplitude = 2;
amplification = 1;
MeasurementAmplification = 1;
u1 = @(t,x,y) amplification* (x.^2 + y.^2 + 1) .* (cos(omega1 .* t) + 2);
u2 = @(t,x,y) amplification* (x.^2 + y.^2 + 1) .* (cos(omega2 .* t) + 2);
u3 = @(t,x,y) u3Amplitude .* u1(t,x,y);

u1grad = @(t,x,y)  amplification.*cat(3,...
    (2.*x).* (cos(omega1 .* t) + 2), ...
    (2.*y).* (cos(omega1 .* t) + 2));
u1F = @(x,y) amplification.* (x.^2 + y.^2 + 1);
u1C = @(x,y) amplification.* (x.^2 + y.^2 + 1) .* 2;
u1Fgrad = @(x,y) amplification.*cat(3,...
    2 .* x, ...
    2 .* y);
u1Cgrad = @(x,y) amplification.*2.*u1Fgrad(x,y);

u1laplace = @(t,x,y) amplification.*4.* (cos(omega1 .* t) + 2);
u1tt = @(t,x,y) amplification.*(x.^2 + y.^2 + 1) .* ((-1).*omega1.^2.*cos(omega1 .* t) );
u1sqtt = @(t,x,y) amplification.^2.*(-2).*omega1.^2.*(x.^2 + y.^2 + 1).^2.*(2.*cos(omega1.*t) + cos(2.*omega1.*t));

u2laplace = @(t,x,y) amplification.* 4.* (cos(omega2 .* t) + 2);
u2tt = @(t,x,y) amplification.*(x.^2 + y.^2 + 1) .* ((-1).*omega2.^2.*cos(omega2 .* t) );
u2sqtt = @(t,x,y) amplification.^2.*(-2).*omega2.^2.*(x.^2 + y.^2 + 1).^2.*(2.*cos(omega2.*t) + cos(2.*omega2.*t));

u3laplace = @(t,x,y) amplification.*u3Amplitude.*4.* (cos(omega3 .* t) + 2);
u3tt = @(t,x,y) amplification.*u3Amplitude.*(x.^2 + y.^2 + 1) .* ((-1).*omega3.^2.*cos(omega3 .* t) );
u3sqtt = @(t,x,y) amplification.^2.*u3Amplitude.^2.*(-2).*omega3.^2.*(x.^2 + y.^2 + 1).^2.*(2.*cos(omega3.*t) + cos(2.*omega3.*t));


% specify our time space cylinder and calculate the triangle mesh in space
% and the mesh in time
timeMeshh = 1/(2*max([f1,f2,f3]));
% time mesh for f1
timeMeshf1 = linspace(0, 1/f1, 1/timeMeshh - 1);
% time mesh for f2
timeMeshf2 = linspace(0, 1/f2, 1/timeMeshh - 1);
% time mesh for f3
timeMeshf3 = linspace(0, 1/f3, 1/timeMeshh - 1);

timeMesh.timeMesh1 = timeMeshf1;
timeMesh.timeMesh2 = timeMeshf2;
timeMesh.timeMesh3 = timeMeshf3;

% our domain
bcenter = [0,0];
brad = 0.2;
domain = [bcenter, brad];

% specify the mesh parameter
meshSize = 0.01;

% compute the triangle mesh
[elements] = initializeMultiLeveLSolver(meshSize, domain);

% prepare FEM matrices a priori
[elements.M_t, elements.tBM,  elements.K, elements.rowK, elements.colK] = prepareFEMMatrices(elements);

% pre-compute some useful things w.r.t. the triangular mesh
fullboundaryIdx = elements.edges(:,1);
interiorIdx = setdiff(1:(size(elements.points,1)), fullboundaryIdx);
elements.interiorIdx = interiorIdx;
elements.boundaryIdx = fullboundaryIdx;
elements.boundaryNormals = 1./sqrt(sum(elements.points(fullboundaryIdx,:).^2,2)).*elements.points(fullboundaryIdx,:); % our center is (0,0), so -> normalisation suffices

% specify the parameters we want to reconstruct
% boundary parameters (not reconstructed)
gamma = 1;

beta = 0;   % this is important check paper for clarification

N = 6; % number of harmonics-1 we will compute
nIter = 6;

% create the space dependent parameters
sourceValueDomain = 2; % B/A of domain

eta_values = [0.01]; % B/A of phantoms
eta_radii = [0.03];
eta_centers = [0.0;0.1];

s_values = [145005]; % B/A of phantoms
s_radii = [0.03];
s_centers = [0.1;-0.1];

b_values = [701]; % B/A of phantoms
b_radii = [0.03];
b_centers = [-0.1;-0.1];

%% so so working
% specify our reference states
f1 = 60;    % Hz
f2 = 78;    % Hz
f3 = f1;    % frequency of third reference state = frequency of first reference state
omega1 = 2*pi*f1;
omega2 = 2*pi*f2;
omega3 = 2*pi*f3;
u3Amplitude = 4;
amplification = 1;
MeasurementAmplification = 1;
u1 = @(t,x,y) amplification* (x.^2 + y.^2 + 1) .* (cos(omega1 .* t) + 2);
u2 = @(t,x,y) amplification* (x.^2 + y.^2 + 1) .* (cos(omega2 .* t) + 2);
u3 = @(t,x,y) u3Amplitude .* u1(t,x,y);

u1grad = @(t,x,y)  amplification.*cat(3,...
    (2.*x).* (cos(omega1 .* t) + 2), ...
    (2.*y).* (cos(omega1 .* t) + 2));
u1F = @(x,y) amplification.* (x.^2 + y.^2 + 1);
u1C = @(x,y) amplification.* (x.^2 + y.^2 + 1) .* 2;
u1Fgrad = @(x,y) amplification.*cat(3,...
    2 .* x, ...
    2 .* y);
u1Cgrad = @(x,y) amplification.*2.*u1Fgrad(x,y);

u1laplace = @(t,x,y) amplification.*4.* (cos(omega1 .* t) + 2);
u1tt = @(t,x,y) amplification.*(x.^2 + y.^2 + 1) .* ((-1).*omega1.^2.*cos(omega1 .* t) );
u1sqtt = @(t,x,y) amplification.^2.*(-2).*omega1.^2.*(x.^2 + y.^2 + 1).^2.*(2.*cos(omega1.*t) + cos(2.*omega1.*t));

u2laplace = @(t,x,y) amplification.* 4.* (cos(omega2 .* t) + 2);
u2tt = @(t,x,y) amplification.*(x.^2 + y.^2 + 1) .* ((-1).*omega2.^2.*cos(omega2 .* t) );
u2sqtt = @(t,x,y) amplification.^2.*(-2).*omega2.^2.*(x.^2 + y.^2 + 1).^2.*(2.*cos(omega2.*t) + cos(2.*omega2.*t));

u3laplace = @(t,x,y) amplification.*u3Amplitude.*4.* (cos(omega3 .* t) + 2);
u3tt = @(t,x,y) amplification.*u3Amplitude.*(x.^2 + y.^2 + 1) .* ((-1).*omega3.^2.*cos(omega3 .* t) );
u3sqtt = @(t,x,y) amplification.^2.*u3Amplitude.^2.*(-2).*omega3.^2.*(x.^2 + y.^2 + 1).^2.*(2.*cos(omega3.*t) + cos(2.*omega3.*t));


% specify our time space cylinder and calculate the triangle mesh in space
% and the mesh in time
timeMeshh = 1/(2*max([f1,f2,f3]));
% time mesh for f1
timeMeshf1 = linspace(0, 1/f1, 1/timeMeshh - 1);
% time mesh for f2
timeMeshf2 = linspace(0, 1/f2, 1/timeMeshh - 1);
% time mesh for f3
timeMeshf3 = linspace(0, 1/f3, 1/timeMeshh - 1);

timeMesh.timeMesh1 = timeMeshf1;
timeMesh.timeMesh2 = timeMeshf2;
timeMesh.timeMesh3 = timeMeshf3;

% our domain
bcenter = [0,0];
brad = 0.2;
domain = [bcenter, brad];

% specify the mesh parameter
meshSize = 0.01;

% compute the triangle mesh
[elements] = initializeMultiLeveLSolver(meshSize, domain);

% prepare FEM matrices a priori
[elements.M_t, elements.tBM,  elements.K, elements.rowK, elements.colK] = prepareFEMMatrices(elements);

% pre-compute some useful things w.r.t. the triangular mesh
fullboundaryIdx = elements.edges(:,1);
interiorIdx = setdiff(1:(size(elements.points,1)), fullboundaryIdx);
elements.interiorIdx = interiorIdx;
elements.boundaryIdx = fullboundaryIdx;
elements.boundaryNormals = 1./sqrt(sum(elements.points(fullboundaryIdx,:).^2,2)).*elements.points(fullboundaryIdx,:); % our center is (0,0), so -> normalisation suffices

% specify the parameters we want to reconstruct
% boundary parameters (not reconstructed)
gamma = 1;

beta = 0;   % this is important check paper for clarification

N = 6; % number of harmonics-1 we will compute
nIter = 6;

% create the space dependent parameters
sourceValueDomain = 2; % B/A of domain

eta_values = [0.01]; % B/A of phantoms
eta_radii = [0.03];
eta_centers = [0.0;0.1];

s_values = [1454]; % B/A of phantoms
s_radii = [0.03];
s_centers = [0.1;-0.1];

b_values = [20.05]; % B/A of phantoms
b_radii = [0.03];
b_centers = [-0.1;-0.1];

%%
% create the space dependent parameters
sourceValueDomain = 2; % B/A of domain

eta_values = [0.001]; % B/A of phantoms
eta_radii = [0.05];
eta_centers = [0.0;0.1];

s_values = [2.5]; % B/A of phantoms
s_radii = [0.05];
s_centers = [0.1;0.1];

b_values = [1.005]; % B/A of phantoms
b_radii = [0.05];
b_centers = [-0.1;-0.1];

eta = constructParameter(elements, eta_centers, eta_radii, eta_values,0);
s = constructParameter(elements, s_centers, s_radii, s_values,1);
b = constructParameter(elements, b_centers, b_radii, b_values,1);