clear all

% specify our reference states
f1 = 25;    % 10 Hz
f2 = 43;    % 20 Hz
f3 = f1;
omega1 = 2*pi*f1;
omega2 = 2*pi*f2;
omega3 = 2*pi*f3;
u3Amplitude = 1.5;
amplification = 0.5;


u1 = @(t,x,y) amplification*(x.^2 + y.^2 + 1) .* (cos(omega1 .* t) + 2);
u2 = @(t,x,y) amplification*(x.^2 + y.^2 + 1) .* (cos(omega2 .* t) + 2);
u3 = @(t,x,y) amplification.*u3Amplitude .* u1(t,x,y);

u1grad = @(t,x,y)  amplification.*cat(3,...
    (2.*x).* (cos(omega1 .* t) + 2), ...
    (2.*y).* (cos(omega1 .* t) + 2));
u1F = @(x,y) amplification.*(x.^2 + y.^2 + 1);
u1C = @(x,y) amplification.*(x.^2 + y.^2 + 1) .* 2;
u1Fgrad = @(x,y) amplification.*cat(3,...
    2 .* x, ...
    2 .* y);
u1Cgrad = @(x,y) amplification.*2.*u1Fgrad(x,y); 

u1laplace = @(t,x,y) amplification.*4.* (cos(omega1 .* t) + 2);
u1tt = @(t,x,y) amplification.*(x.^2 + y.^2 + 1) .* ((-1).*omega1.^2.*cos(omega1 .* t) );
u1sqtt = @(t,x,y) amplification.^2.*(-2).*omega1.^2.*(x.^2 + y.^2 + 1).^2.*(2.*cos(omega1.*t) + cos(2.*omega1.*t));

u2laplace = @(t,x,y) amplification.* 4.* (cos(omega2 .* t) + 2);
u2tt = @(t,x,y) amplification.*(x.^2 + y.^2 + 1) .* ((-1).*omega2.^2.*cos(omega2 .* t) );
u2sqtt = @(t,x,y) amplification.^2.*(-2).*omega2.^2.*(x.^2 + y.^2 + 1).^2.*(2.*cos(omega2.*t)+ cos(2.*omega2.*t));

u3laplace = @(t,x,y) amplification.*u3Amplitude.*4.* (cos(omega3 .* t) + 2);
u3tt = @(t,x,y) amplification.*u3Amplitude.*(x.^2 + y.^2 + 1) .* ((-1).*omega3.^2.*cos(omega3 .* t) );
u3sqtt = @(t,x,y) amplification.^2.*u3Amplitude.^2.*(-2).*omega3.^2.*(x.^2 + y.^2 + 1).^2.*(2.*cos(omega3.*t)+ cos(2.*omega3.*t));


% specify our time space cylinder and calculate the triangle mesh in space
% and the mesh in time
timeMeshh = 0.01;
% time domain (lowest frequency determines the duration)
timeMesh = linspace(0, 1/f1, 1/timeMeshh - 1);
% recompute time difference
timeMeshh = timeMesh(2) - timeMesh(1); % careful this is the time diff!

% our domain
bcenter = [0,0];
brad = 0.2;
domain = [bcenter, brad];

% specify the mesh parameter
meshSize = 0.01;

% compute the triangle mesh
[elements] = initializeMultiLeveLSolver(meshSize, domain);

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

% define a phantom in our domain with different speed of sound, diffusivity
% and nonlinearity parameter
diffusivity = 0.5;
values = [3, 3]; % B/A of phantoms
radii = [0.05,0.05];
diffusivityPhantoms = [0.49,0.49]; % this allows to adjust the diffusivity for the phantoms
centers = [0,0; 0.1,-0.1];

massDensity = 1000; %kg/m^3

speed_of_sound = 2;

N = 6; % number of harmonics-1 we will compute
nIter = 6;

% create the space dependent parameters
sourceValueDomain = 2; % B/A of domain

eta = constructNonlinearityDivB(elements, massDensity, speed_of_sound, diffusivity, diffusivityPhantoms, centers, radii, values, sourceValueDomain, true); %nonlinearity scaled by 1/b

s = constructSquaredSpeedOfSoundDivB(elements, speed_of_sound, diffusivity, diffusivityPhantoms, centers, radii); % speed of sound scaled by 1/b

b = constructReciprocalDiffusivity(elements, diffusivity, diffusivityPhantoms, centers, radii);

% the complex wavenumber, here we compute the square wave number
kappasq = constructKappaReparameterized(elements, s, b, [omega1 omega2 omega3], N); % compute all the complex wave numbers needed

% since u^0_j  does not need to solve the PDE,...
% construct the boundary excitations (same source just with different
% frequency

%% prepare the source(s)
sourceEdge = 1; % we impose the source on the boundary (negative quadrant)
%boundaryPointsSourceIdx = elements.edges((elements.edges(:,3) == sourceEdge),1);

boundaryPointsSourceIdx = elements.boundaryIdx;

boundaryPointsSource = elements.points(boundaryPointsSourceIdx,:);

% compute the normals
boundaryPointsSourceNormals = 1./sqrt(sum(elements.points(boundaryPointsSourceIdx,:).^2,2)).*elements.points(boundaryPointsSourceIdx,:); % our center is (0,0), so -> normalisation is suffices

% frequency part (robin boundary)
sourceFrequency = zeros(size(elements.points,1),1);
sourceFrequency(boundaryPointsSourceIdx) = (gamma.*u1F(boundaryPointsSource(:,1),boundaryPointsSource(:,2)) + dot(squeeze(u1Fgrad(boundaryPointsSource(:,1),boundaryPointsSource(:,2))).',boundaryPointsSourceNormals.').');

% constant part (robin boundary)
sourceConstant = zeros(size(elements.points,1),1);
sourceConstant(boundaryPointsSourceIdx) = gamma.*u1C(boundaryPointsSource(:,1), boundaryPointsSource(:,2)) + dot(squeeze(u1Cgrad(boundaryPointsSource(:,1),boundaryPointsSource(:,2))).',boundaryPointsSourceNormals.').';

%%
excitations = zeros(size(elements.points,1), N, 3);
excitations(:,1,1) = amplification.*sourceConstant;
excitations(:,2,1) = amplification.*sourceFrequency;
excitations(:,1,2) = amplification.*sourceConstant;
excitations(:,2,2) = amplification.*sourceFrequency;
excitations(:,1,3) = amplification.*u3Amplitude.*sourceConstant;
excitations(:,2,3) = amplification.*u3Amplitude.*sourceFrequency;
%%
[cN, U1, F1] = solveWesterveltMultiLevelBoundaryExcitation(elements, omega1, beta, gamma, squeeze(kappasq(:,:,1)), squeeze(excitations(:,:,1)), eta, b, nIter, N, 10^(-12));
[cN, U2, F2] = solveWesterveltMultiLevelBoundaryExcitation(elements, omega2, beta, gamma, squeeze(kappasq(:,:,2)), squeeze(excitations(:,:,2)), eta, b, nIter, N, 10^(-12));
[cN, U3, F3] = solveWesterveltMultiLevelBoundaryExcitation(elements, omega3, beta, gamma, squeeze(kappasq(:,:,3)), squeeze(excitations(:,:,3)), eta, b, nIter, N, 10^(-12));
%%
% compute the solutions on the time - space mesh
u1s = calcSolution(timeMesh, squeeze(U1(N,:,:)), omega1);
u2s = calcSolution(timeMesh, squeeze(U2(N,:,:)), omega2);
u3s = calcSolution(timeMesh, squeeze(U3(N,:,:)), omega3);

%%
% define \Sigma our measurement manifold/discrete points
% the triangulation already defines the the edges of our doimain 1:4
% (circle sectors)

measurementEdge = 3; % positive quadrant edge

% fetch the boundary points
elements.measurementPointsIdx = elements.edges((elements.edges(:,3) == measurementEdge),1);

%elements.measurementPointsIdx = measurementPointsIdx;

% do the measurement on the whole boundary
elements.measurementPointsIdx = elements.boundaryIdx;

measurement_u1 = u1s(:,elements.measurementPointsIdx);
measurement_u2 = u2s(:,elements.measurementPointsIdx);
measurement_u3 = u3s(:,elements.measurementPointsIdx);

measurement_u1_harmonics = zeros(size(squeeze(U1(N,:,:))));

measurement_u2_harmonics = zeros(size(squeeze(U2(N,:,:))));

measurement_u3_harmonics = zeros(size(squeeze(U3(N,:,:))));


for j=1:N
    measurement_u1_harmonics(j,elements.measurementPointsIdx) = squeeze(U1(N,j,elements.measurementPointsIdx)).';

    measurement_u2_harmonics(j,elements.measurementPointsIdx) = squeeze(U2(N,j,elements.measurementPointsIdx)).';

    measurement_u3_harmonics(j,elements.measurementPointsIdx) = squeeze(U3(N,j,elements.measurementPointsIdx)).';
end


%%
% pre compoute the gradient operator
p = elements.points;
t = elements.tri(:,1:3);

x = p(:,1);
y = p(:,2);

nt = size(t,1);
np = size(p,1);

% triangle vertex indices
i1 = t(:,1); i2 = t(:,2); i3 = t(:,3);

% coordinates
x1 = x(i1); y1 = y(i1);
x2 = x(i2); y2 = y(i2);
x3 = x(i3); y3 = y(i3);

% triangle areas
area = 0.5*((x2-x1).*(y3-y1)-(x3-x1).*(y2-y1));

% gradients of basis functions
b1 = (y2-y3)./(2*area);
b2 = (y3-y1)./(2*area);
b3 = (y1-y2)./(2*area);

c1 = (x3-x2)./(2*area);
c2 = (x1-x3)./(2*area);
c3 = (x2-x1)./(2*area);

% assemble sparse gradient matrices
rows = repmat((1:nt)',1,3);

Gx_elem = sparse(rows(:), t(:), [b1;b2;b3], nt, np);
Gy_elem = sparse(rows(:), t(:), [c1;c2;c3], nt, np);

% average to nodes (mass lumping style)
Mmap = sparse(repmat((1:nt)',3,1), t(:), 1, nt, np);
weights = Mmap' * abs(area);

Gx = (Mmap' * (abs(area).*Gx_elem)) ./ weights;
Gy = (Mmap' * (abs(area).*Gy_elem)) ./ weights;

% sanity check for whether \gamma*u + \nabla u \cdot n = g on the boundary

%testu1boundary = zeros(1,size(u2s,2));
testu1boundary(1,boundaryPointsSourceIdx) = gamma.*u1(0, boundaryPointsSource(:,1), boundaryPointsSource(:,2)) + dot(squeeze(u1grad(0,boundaryPointsSource(:,1),boundaryPointsSource(:,2))).', boundaryPointsSourceNormals.').';

%% for the inverse problem we construct u0_1, u0_2, u0_3 (in this case as solution of our approx. scheme)
% at x0 = (s^0, b^0, 0), where s^0, b^0 are space constant functions
% for finer triangular meshes (high accuracy) these can be precomputed and
% stored to speed up computation
s0 = min(s).*ones(size(s));
b0 = min(b).*ones(size(b));
eta0 = zeros(size(eta));
kappasq0 = constructKappaReparameterized(elements, s0, b0, [omega1 omega2 omega3], N); % compute all the complex wave numbers needed
excitationsReferenceState = excitations;
% for i=1:size(elements.points,1)
%     u1Sampled(i) = u1F(elements.points(i,1), elements.points(i,2));
% end
% the boundary excitation we take from our exemplary function
[~, U0_1, U0_F1] = solveWesterveltMultiLevelBoundaryExcitation(elements, omega1, beta, gamma, squeeze(kappasq0(:,:,1)), squeeze(excitationsReferenceState(:,:,1)), eta0, b0, nIter, N, 10^(-12));
[~, U0_2, U0_F2] = solveWesterveltMultiLevelBoundaryExcitation(elements, omega2, beta, gamma, squeeze(kappasq0(:,:,2)), squeeze(excitationsReferenceState(:,:,2)), eta0, b0, nIter, N, 10^(-12));
[~, U0_3, U0_F3] = solveWesterveltMultiLevelBoundaryExcitation(elements, omega3, beta, gamma, squeeze(kappasq0(:,:,3)), squeeze(excitationsReferenceState(:,:,3)), eta0, b0, nIter, N, 10^(-12));

% there reference state in time x space
u0_1 = calcSolution(timeMesh, squeeze(U0_1(N,:,:)), omega1);
u0_2 = calcSolution(timeMesh, squeeze(U0_2(N,:,:)), omega2);
u0_3 = calcSolution(timeMesh, squeeze(U0_3(N,:,:)), omega3);

% just a quick sanity check --> the solutions MUST fulfill |\Delta u^0|
% \geq c > 0 and |u^0| \geq c > 0 a.e. in our time space cylinder, if we
% are calculating the effective increments, if not, then bäh
if min(min(u0_1)) <= 0 || min(min(u0_2)) <= 0 || min(min(u0_3)) <= 0
    warning('One of the reference states is not bounded from below (ignore this if we do not computing effective increments).');
end
% % check the laplacian
% if min(min(M\L*u0_1.'))<= 0 || min(min(M\L*u0_2.')) <= 0 || min(min(M\L*u0_3.')) <= 0
%     error('Laplacian of one of the reference states is not bounded from below.');
% end

% another sanity check --> check whether the boundary condition is
% fulfilled
u0_1_rb = calcRobinBoundary(elements, u0_1, gamma, Gx, Gy);
u0_2_rb = calcRobinBoundary(elements, u0_2, gamma, Gx, Gy);
u0_3_rb = calcRobinBoundary(elements, u0_3, gamma, Gx, Gy);

for i=1:size(timeMesh,2)

    u1sampled(i,:) = u1(timeMesh(i), elements.points(:,1), elements.points(:,2));
    u2sampled(i,:) = u2(timeMesh(i), elements.points(:,1), elements.points(:,2));
    u3sampled(i,:) = u3(timeMesh(i), elements.points(:,1), elements.points(:,2));

    u1LaplacianSampled(i,:) = u1laplace(timeMesh(i), elements.points(:,1), elements.points(:,2));
    u1ttSampled(i,:) = u1tt(timeMesh(i), elements.points(:,1), elements.points(:,2));
    u1sqttSampled(i,:) = u1sqtt(timeMesh(i), elements.points(:,1), elements.points(:,2));

    u2LaplacianSampled(i,:) = u2laplace(timeMesh(i), elements.points(:,1), elements.points(:,2));
    u2ttSampled(i,:) = u2tt(timeMesh(i), elements.points(:,1), elements.points(:,2));
    u2sqttSampled(i,:) = u2sqtt(timeMesh(i), elements.points(:,1), elements.points(:,2));


    u3LaplacianSampled(i,:) = u3laplace(timeMesh(i), elements.points(:,1), elements.points(:,2));
    u3ttSampled(i,:) = u3tt(timeMesh(i), elements.points(:,1), elements.points(:,2));
    u3sqttSampled(i,:) = u3sqtt(timeMesh(i), elements.points(:,1), elements.points(:,2));

end

u1b = calcRobinBoundary(elements, u1sampled, gamma, Gx, Gy);
u2b = calcRobinBoundary(elements, u2sampled, gamma, Gx, Gy);
u3b = calcRobinBoundary(elements, u3sampled, gamma, Gx, Gy);

u1s_b = calcRobinBoundary(elements, u1s, gamma, Gx, Gy);
u2s_b = calcRobinBoundary(elements, u2s, gamma, Gx, Gy);
u3s_b = calcRobinBoundary(elements, u3s, gamma, Gx, Gy);


% if these values are big, or do not decrease with a finer triangular mesh
% something is odd
l2boundaryValErroru0_1 = sqrt(sum(sum(abs(u1b(:,elements.boundaryIdx) - u0_1_rb(:,elements.boundaryIdx)).^2,2),1));
linfboundaryValErroru0_1 = sqrt(max(max(abs(u1b(:,elements.boundaryIdx) - u0_1_rb(:,elements.boundaryIdx)))));

l2boundaryValErroru0_2 = sqrt(sum(sum(abs(u2b(:,elements.boundaryIdx) - u0_2_rb(:,elements.boundaryIdx)).^2,2),1));
linfboundaryValErroru0_2 = sqrt(max(max(abs(u2b(:,elements.boundaryIdx) - u0_2_rb(:,elements.boundaryIdx)))));

l2boundaryValErroru0_3 = sqrt(sum(sum(abs(u3b(:,elements.boundaryIdx) - u0_3_rb(:,elements.boundaryIdx)).^2,2),1));
linfboundaryValErroru0_3 = sqrt(max(max(abs(u3b(:,elements.boundaryIdx) - u0_3_rb(:,elements.boundaryIdx)))));

% these distances have to be small!
l2boundaryValErroru0_r1 = sqrt(sum(sum(abs(u1b(:,elements.boundaryIdx) - u1s_b(:,elements.boundaryIdx)).^2,2),1));
linfboundaryValErroru0_r1 = sqrt(max(max(abs(u1b(:,elements.boundaryIdx) - u1s_b(:,elements.boundaryIdx)))));

l2boundaryValErroru0_r2 = sqrt(sum(sum(abs(u2b(:,elements.boundaryIdx) - u2s_b(:,elements.boundaryIdx)).^2,2),1));
linfboundaryValErroru0_r2 = sqrt(max(max(abs(u2b(:,elements.boundaryIdx) - u2s_b(:,elements.boundaryIdx)))));

l2boundaryValErroru0_r3 = sqrt(sum(sum(abs(u3b(:,elements.boundaryIdx) - u3s_b(:,elements.boundaryIdx)).^2,2),1));
linfboundaryValErroru0_r3 = sqrt(max(max(abs(u3b(:,elements.boundaryIdx) - u3s_b(:,elements.boundaryIdx)))));

% check if x0 is really close enough to the solution (in L^2(\Omega))
[~, s0error] = integrate_fun_trimesh(elements.opoints, elements.otri, ((s0 - s).^2).');
[~, b0error] = integrate_fun_trimesh(elements.opoints, elements.otri, ((b0 - b).^2).');
[~, eta0error] = integrate_fun_trimesh(elements.opoints, elements.otri, ((eta).^2).'); % eta0 = 0, always!

errorToSol = sqrt( s0error + b0error + eta0error);

if errorToSol > 1
    warning("Initial x0 is not near enough to the solution.");
end

%% if the error in the parameters is small we should have F(x_0) + DF(x_0)(x-x_0) \approx F(x)

% first the obvious one

%F(x) + DF(x_0)(x) = F(x)
dx.ds = 0;
dx.db = 0;
dx.deta = 0;
dx.excitation = zeros(size(elements.points,1), N);

x0.refState_1.s0 = s0;
x0.refState_1.b0 = b0;
x0.refState_1.eta0 = eta0;
x0.refState_1.u0 = squeeze(U0_1(N,:,:));
x0.refState_1.F = U0_F1;
x0.refState_1.kappa0 = constructKappaReparameterized(elements, x0.refState_1.s0, x0.refState_1.b0, omega1, N);

[~, DU, ~] = solveLinearisedWesterveltMultiLevelBoundaryExcitation(elements, omega1, beta, gamma, x0.refState_1, dx, nIter, N, true);

du = calcSolution(timeMesh, squeeze(DU(N,:,:)), omega1);

err = u1s + du - u1s;
if max(max(abs(err))) > 10e-15
    error("DF(x_0)(0) is not zero");
end

% now with a small pertubation
perturbationCoeff = 0.9;
perturbed_s = s.*perturbationCoeff;
perturbed_b = b.*perturbationCoeff;
perturbed_eta = 0;

kappaPerturbed = constructKappaReparameterized(elements, perturbed_s, perturbed_b, omega1, N);

[~, UPert, FPert] = solveWesterveltMultiLevelBoundaryExcitation(elements, omega1, beta, gamma, kappaPerturbed, squeeze(excitationsReferenceState(:,:,1)), perturbed_eta, perturbed_b, nIter, N, 10^(-12));

dx.ds = s - perturbed_s;
dx.db = b - perturbed_b;
dx.deta = eta  - perturbed_eta;
dx.excitation = zeros(size(elements.points,1), N);

xs = x0;
xs.refState_1.u0 = squeeze(UPert(N,:,:));
xs.refState_1.F = FPert;
xs.refState_1.s0 = perturbed_s;
xs.refState_1.b0 = perturbed_b;
xs.refState_1.eta0 = perturbed_eta;

[~, DU, ~] = solveLinearisedWesterveltMultiLevelBoundaryExcitation(elements, omega1, beta, gamma, xs.refState_1, dx, nIter, N, true);
du = calcSolution(timeMesh, squeeze(DU(N,:,:)), omega1);
[~, ldx.ds] = integrate_fun_trimesh(elements.opoints, elements.otri, abs(dx.ds).^2.');
[~, ldx.db] = integrate_fun_trimesh(elements.opoints, elements.otri, abs(dx.db).^2.');
[~, ldx.deta] = integrate_fun_trimesh(elements.opoints, elements.otri, abs(dx.deta).^2.');

nm = sqrt(ldx.ds + ldx.db + ldx.deta);
% for very small pertubations we should have F(x) \approx F(x0)
diffU = (UPert(N,:,:) - DU(N,:,:) - U1(N,:,:));
diffU_time = calcSolution(timeMesh, squeeze(diffU), omega1);
for i = 1:size(timeMesh,2)
    [~, a(i)] = integrate_fun_trimesh(elements.opoints, elements.otri, abs(diffU_time(1,:)).^2);
end
dn = trapz(timeMesh,a);

%% check the deviation in ds only 
perturbationCoeff = 0.7;
perturbed_s = s.*perturbationCoeff;
perturbed_b = b0;
perturbed_eta = eta0;

kappaPerturbed = constructKappaReparameterized(elements, perturbed_s, perturbed_b, omega1, N);

[~, UPert, FPert] = solveWesterveltMultiLevelBoundaryExcitation(elements, omega1, beta, gamma, kappaPerturbed, squeeze(excitationsReferenceState(:,:,1)), perturbed_eta, perturbed_b, nIter, N, 10^(-12));

dx.ds = s - perturbed_s;
dx.db = b - perturbed_b;
dx.deta = eta  - perturbed_eta;
dx.excitation = zeros(size(elements.points,1), N);

xs = x0;
xs.refState_1.u0 = squeeze(UPert(N,:,:));
xs.refState_1.F = FPert;
xs.refState_1.s0 = perturbed_s;
xs.refState_1.b0 = perturbed_b;
xs.refState_1.eta0 = perturbed_eta;

[~, DU, ~] = solveLinearisedWesterveltMultiLevelBoundaryExcitation(elements, omega1, beta, gamma, xs.refState_1, dx, nIter, N, true);
du = calcSolution(timeMesh, squeeze(DU(N,:,:)), omega1);
[~, ldx.ds] = integrate_fun_trimesh(elements.opoints, elements.otri, abs(dx.ds).^2.');
[~, ldx.db] = integrate_fun_trimesh(elements.opoints, elements.otri, abs(dx.db).^2.');
[~, ldx.deta] = integrate_fun_trimesh(elements.opoints, elements.otri, abs(dx.deta).^2.');

nm = sqrt(ldx.ds + ldx.db + ldx.deta);
% for very small pertubations we should have F(x) \approx F(x0)
diffU = (UPert(N,:,:) - DU(N,:,:) - U1(N,:,:));
diffU_time = calcSolution(timeMesh, squeeze(diffU), omega1);
for i = 1:size(timeMesh,2)
    [~, a(i)] = integrate_fun_trimesh(elements.opoints, elements.otri, abs(diffU_time(1,:)).^2);
end
dn = trapz(timeMesh,a);
harmonics = 0:(N-1);
lapu0 = squeeze(kappasq0(:,:,1)).'.*squeeze(U0_1(N,:,:));
calcLap0sol = calcSolution(timeMesh, -omega1.^2.*harmonics.'.^2 .* lapu0, omega1);
% calc also the adjoint state
residual_1 = zeros(N,size(elements.points,1));
residual_1(:,elements.measurementPointsIdx) = squeeze(DU(N,:,elements.measurementPointsIdx));
[~, Uadj_1, ~] = solveAdjointLinearisedWesterveltMultiLevelBoundaryExcitation(elements, omega1, beta, gamma, xs.refState_1, residual_1, nIter, N);
uadj = calcSolution(timeMesh, squeeze(Uadj_1(N,:,:)), omega1);
ads = trapz(timeMesh, uadj.*(calcLap0sol));

%% theory tells us that we do not need that u0 is a solution of our PDE
% prepare the harmonics of our reference states
u0Csampled = u1C(elements.points(:,1), elements.points(:,2));
u0Fsampled = u1F(elements.points(:,1), elements.points(:,2));
u0sampled = zeros(size(squeeze(U0_1(cN,:,:))));
laplaceu0 = zeros(size(squeeze(U0_1(cN,:,:))));
u0sampled(1,:) = u0Csampled;
u0sampled(2,:) = u0Fsampled;
laplaceu0(1,:) = 8;
laplaceu0(2,:) = 4;
% check whether our fourier transform is correct
u0sampledrecon = calcSolution(timeMesh, u0sampled, omega1);

% another sanity check
if norm(norm(abs(u0sampledrecon - u1sampled),2),2) > 10e-8
    error('Fourier coefficients of reference state 1 do not match.');
end

% check also the higher amplitude reference state
u0_3sampledrecon = calcSolution(timeMesh, amplification.*u3Amplitude*u0sampled, omega1);

% another sanity check
if norm(norm(abs(u0_3sampledrecon - u3sampled),2),2) > 10e-8
    error('Fourier coefficients of reference state 3 do not match.');
end

%%
useSolutionAsLinPoint = false;


if useSolutionAsLinPoint == true

    x0.refState_1.u0 = squeeze(U0_1(N,:,:));
    x0.refState_1.F = U0_F1;
    x0.refState_1.kappa0 = constructKappaReparameterized(elements, s0, b0, omega1, N);

    x0.refState_2.u0 = squeeze(U0_2(N,:,:));
    x0.refState_2.F = U0_F2;
    x0.refState_2.kappa0 = constructKappaReparameterized(elements, s0, b0, omega2, N);

    x0.refState_3.u0 = squeeze(U0_3(N,:,:));
    x0.refState_3.F = U0_F3;
    x0.refState_3.kappa0 = constructKappaReparameterized(elements, s0, b0, omega3, N);

    for j = 0:(N-1)
        u1laplacef(j+1,:) = -x0.refState_1.kappa0(:,j+1).'.*j^2.*x0.refState_1.u0(j+1,:) - x0.refState_1.F(j+1,:);
        u1ttf(j+1,:)  = -j^2.*omega1.^2.*x0.refState_1.u0(j+1,:);

        u2laplacef(j+1,:) = -x0.refState_2.kappa0(:,j+1).'.*j^2.*x0.refState_2.u0(j+1,:) - x0.refState_2.F(j+1,:);
        u2ttf(j+1,:)  = -j^2.*omega2.^2.*x0.refState_2.u0(j+1,:);

        u3laplacef(j+1,:) = -x0.refState_3.kappa0(:,j+1).'.*j^2.*x0.refState_3.u0(j+1,:) - x0.refState_3.F(j+1,:);
        u3ttf(j+1,:)  = -j^2.*omega3.^2.*x0.refState_3.u0(j+1,:);

        % computing u^2_{tt} is a bit more tricky
        p_m = zeros(1,size(elements.points,1));
        p_m2 = zeros(1,size(elements.points,1));
        p_m3 = zeros(1,size(elements.points,1));

        for l = 0:j
            p_m = p_m + ...
                squeeze(x0.refState_1.u0(l+1,:)) .* ...
                squeeze(x0.refState_1.u0((j-l)+1,:));

            p_m2 = p_m2 + ...
                squeeze(x0.refState_2.u0(l+1,:)) .* ...
                squeeze(x0.refState_2.u0((j-l)+1,:));

            p_m3 = p_m3 + ...
                squeeze(x0.refState_3.u0(l+1,:)) .* ...
                squeeze(x0.refState_3.u0((j-l)+1,:));
        end

        % ---------- Second sum ----------
        % 2 * sum_{r=0}^{N-1-j} conj(u_r) * u_{r+j}
        for r = j:2:(2*(N-1) - j)
            minusidx = (r-j)/2;
            plusidx = (r+j)/2;
            p_m = p_m + 2 * ...
                conj(squeeze(x0.refState_1.u0(minusidx+1,:))) .* ...
                squeeze(x0.refState_1.u0(plusidx+1,:));

            p_m2 = p_m2 + 2 * ...
                conj(squeeze(x0.refState_2.u0(minusidx+1,:))) .* ...
                squeeze(x0.refState_2.u0(plusidx+1,:));
            p_m3 = p_m3 + 2 * ...
                conj(squeeze(x0.refState_3.u0(minusidx+1,:))) .* ...
                squeeze(x0.refState_3.u0(plusidx+1,:));
        end
        u1sqttf(j+1,:) = -j.^2.*omega1^2.*p_m;
        u2sqttf(j+1,:) = -j.^2.*omega2^2.*p_m2;
        u3sqttf(j+1,:) = -j.^2.*omega3^2.*p_m3;
    end
    u1tt = calcSolution(timeMesh, u1ttf, omega1);
    u1lap = calcSolution(timeMesh, u1laplacef,omega1);
    u1sqtt = calcSolution(timeMesh, u1sqttf,omega1);

    u2tt = calcSolution(timeMesh, u2ttf, omega2);
    u2lap = calcSolution(timeMesh, u2laplacef,omega2);
    u2sqtt = calcSolution(timeMesh, u2sqttf,omega2);

    u3tt = calcSolution(timeMesh, u3ttf, omega3);
    u3lap = calcSolution(timeMesh, u3laplacef,omega3);
    u3sqtt = calcSolution(timeMesh, u3sqttf,omega3);
    
    referenceStates.u1LaplacianSampled = u1lap;
    referenceStates.u1ttSampled = u1tt;
    referenceStates.u1sqttSampled = u1sqtt;

    referenceStates.u2LaplacianSampled = u2lap;
    referenceStates.u2ttSampled = u2tt;
    referenceStates.u2sqttSampled = u2sqtt;

    referenceStates.u3LaplacianSampled = u3lap;
    referenceStates.u3ttSampled = u3tt;
    referenceStates.u3sqttSampled = u3sqtt;

else
    kappasq0 = constructKappaReparameterized(elements, s0, b0, [omega1 omega2 omega3], N); % compute all the complex wave numbers needed

    x0.refState_1.u0 = u0sampled;
    x0.refState_1.laplaceu0 = laplaceu0;
    x0.refState_1.kappa0 = squeeze(kappasq0(:,:,1));

    x0.refState_2.u0 = u0sampled;
    x0.refState_2.laplaceu0 = laplaceu0;
    x0.refState_2.kappa0 = squeeze(kappasq0(:,:,2));

    x0.refState_3.u0 = amplification.*u3Amplitude*u0sampled;
    x0.refState_3.laplaceu0 = amplification.*u3Amplitude*laplaceu0;
    x0.refState_3.kappa0 = squeeze(kappasq0(:,:,3));

    referenceStates.u1LaplacianSampled = u1LaplacianSampled;
    referenceStates.u1ttSampled = u1ttSampled;
    referenceStates.u1sqttSampled = u1sqttSampled;

    referenceStates.u2LaplacianSampled = u2LaplacianSampled;
    referenceStates.u2ttSampled = u2ttSampled;
    referenceStates.u2sqttSampled = u2sqttSampled;

    referenceStates.u3LaplacianSampled = u3LaplacianSampled;
    referenceStates.u3ttSampled = u3ttSampled;
    referenceStates.u3sqttSampled = u3sqttSampled;

end

x0.refState_1.s0 = s0;
x0.refState_1.b0 = b0;
x0.refState_1.eta0 = eta0;

x0.refState_2.s0 = s0;
x0.refState_2.b0 = b0;
x0.refState_2.eta0 = eta0;

x0.refState_3.s0 = s0;
x0.refState_3.b0 = b0;
x0.refState_3.eta0 = eta0;

xsol = frozenNewtonMethod(elements, timeMesh, x0, referenceStates, beta, gamma, measurement_u1_harmonics, measurement_u2_harmonics, measurement_u3_harmonics, omega1, omega2, omega3, excitations, useSolutionAsLinPoint, nIter, N, 1000, 10e-10, 10e-14);
%%
point = [0.0;0.05];

[v,idx] = min(sum((elements.points - point(:)').^2,2));

node = [elements.points(idx,1);elements.points(idx,2)];

T = 1/f2;
% sampling frequency in time
Fs = 1/T * 2 * (N);
omega = omega1;
U = squeeze(U1(N,:,:));

Ns = 2000;
pC = zeros(1,Ns);
for m=0:(N-1)
    pC = pC + U(m+1,idx) .* exp(1i.*m.*omega.*(0:(Ns-1))*1/Fs);
end

MaxBins = 5;
P0 = max(max(u2s(:,:)));

window = hanning(Ns);
freq =  Fs/Ns*(0:(Ns/2));
y = abs(fft(window'.*real(pC)))/Ns;
y1 = y(1:Ns/2+1);
y1(2:end-1) = 2*y1(2:end-1);
shiftedTF = fftshift(fft(real(pC)))/Ns;
TFdB = 10*log10(y1/P0);
fscaling = 10^3;
M = min(MaxBins, Ns/2 + 1);
xaxis = freq./fscaling;
figure, plot(xaxis, TFdB(1:Ns/2+1))
xlabel("Frequency [kHz]")
ylabel("P/P0 [dB]")
