clear all

% specify our reference states
f1 = 10;    % 10 Hz
f2 = 20;    % 20 Hz             
omega1 = 2*pi*f1;     
omega2 = 2*pi*f2;     

u1 = @(t,x,y) (x.^2 + y.^2 + 1) .* (cos(omega1 .* t) + 2);
u2 = @(t,x,y) (x.^2 + y.^2 + 1) .* (cos(omega2 .* t) + 2);
u3 = @(t,x,y) 2 .* u1(t,x,y);

u1grad = @(t,x,y)  cat(3,...
    (2.*x).* (cos(omega1 .* t) + 2), ...
    (2.*y).* (cos(omega1 .* t) + 2));
u1F = @(x,y) 1./2.*(x.^2 + y.^2 + 1); 
u1C = @(x,y) (x.^2 + y.^2 + 1) .* 2;
u1Fgrad = @(x,y) cat(3,...
    2 .* x, ...
    2 .* y);
u1Cgrad = @(x,y) 2.*u1Fgrad(x,y);

% specify our time space cylinder and calculate the triangle mesh in space
% and the mesh in time
timeMeshh = 0.01;
% time domain (lowest frequency determines the duration)
timeMesh = linspace(0,1/f1,1/timeMeshh - 1);
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
diffusivity = 0.005;
values = [5]; % B/A of phantoms
radii = [0.05];
diffusivityPhantoms = [0.05]; % this allows to adjust the diffusivity for the phantoms
centers = [0; 0];

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
kappasq = constructKappaReparameterized(elements, s, b, [omega1 omega2 omega1], N); % compute all the complex wave numbers needed

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
sourceFrequency(boundaryPointsSourceIdx) = 1/2.*(gamma.*u1F(boundaryPointsSource(:,1),boundaryPointsSource(:,2)) + dot(squeeze(u1Fgrad(boundaryPointsSource(:,1),boundaryPointsSource(:,2))).',boundaryPointsSourceNormals.').');

% constant part (robin boundary)
sourceConstant = zeros(size(elements.points,1),1);
sourceConstant(boundaryPointsSourceIdx) = gamma.*u1C(boundaryPointsSource(:,1), boundaryPointsSource(:,2)) + dot(squeeze(u1Cgrad(boundaryPointsSource(:,1),boundaryPointsSource(:,2))).',boundaryPointsSourceNormals.').';

%%
excitations = zeros(size(elements.points,1), N, 3);
excitations(:,1,1) = sourceConstant;
excitations(:,2,1) = sourceFrequency;
excitations(:,1,2) = sourceConstant;
excitations(:,2,2) = sourceFrequency;
excitations(:,1,3) = 2.*sourceConstant;
excitations(:,2,3) = 2.*sourceFrequency;
%%
[cN, U1, F1] = solveWesterveltMultiLevelBoundaryExcitation(elements, omega1, beta, gamma, squeeze(kappasq(:,:,1)), squeeze(excitations(:,:,1)), eta, b, nIter, N, 10^(-12));
[cN, U2, F2] = solveWesterveltMultiLevelBoundaryExcitation(elements, omega2, beta, gamma, squeeze(kappasq(:,:,2)), squeeze(excitations(:,:,2)), eta, b, nIter, N, 10^(-12));
[cN, U3, F3] = solveWesterveltMultiLevelBoundaryExcitation(elements, omega1, beta, gamma, squeeze(kappasq(:,:,3)), squeeze(excitations(:,:,3)), eta, b, nIter, N, 10^(-12));
%%
% compute the solutions on the time - space mesh

u1s = calcSolution(timeMesh, squeeze(U1(N,:,:)), omega1);
u2s = calcSolution(timeMesh, squeeze(U2(N,:,:)), omega2);
u3s = calcSolution(timeMesh, squeeze(U3(N,:,:)), omega1);


%%
% define \Sigma our measurement manifold/discrete points
% the triangulation already defines the the edges of our doimain 1:4
% (circle sectors)

measurementEdge = 3; % positive quadrant edge

% fetch the boundary points
measurementPointsIdx = elements.edges((elements.edges(:,3) == measurementEdge),1);
elements.measurementPointsIdx = measurementPointsIdx;
measurement_u1 = u1s(:,measurementPointsIdx);
measurement_u2 = u2s(:,measurementPointsIdx);
measurement_u3 = u3s(:,measurementPointsIdx);

measurement_u1_harmonics = zeros(size(squeeze(U1(N,:,:))));

measurement_u2_harmonics = zeros(size(squeeze(U2(N,:,:))));

measurement_u3_harmonics = zeros(size(squeeze(U3(N,:,:))));


for j=1:N
    measurement_u1_harmonics(j,measurementPointsIdx) = squeeze(U1(N,j,measurementPointsIdx)).';

    measurement_u2_harmonics(j,measurementPointsIdx) = squeeze(U2(N,j,measurementPointsIdx)).';

measurement_u3_harmonics(j,measurementPointsIdx) = squeeze(U3(N,j,measurementPointsIdx)).';
end


%%
% There is a bunch of things we can prepare beforehand computation
L = cotmatrix(elements.points, elements.tri); % laplacian matrix (space) [stiffnes matrix]
M = massmatrix(elements.points, elements.tri); % mass matrix (space)

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

%testu1boundary = zeros(1,size(u2s,2));
% testu1boundary(1,boundaryPointsSourceIdx) = gamma.*u1(0, boundaryPointsSource(:,1), boundaryPointsSource(:,2)) + dot(squeeze(u1grad(0,boundaryPointsSource(:,1),boundaryPointsSource(:,2))).', boundaryPointsSourceNormals.').';
% this is the check for the scaled version:
%figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2), real(-squeeze(kappasq(:,2,1)).*squeeze(U1(N,2,:)) - M\L * squeeze(U1(N,2,:))), 'facecolor', 'interp'); shading interp;
% this is the check for the unscaled version:
%figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2), real(-squeeze(-b.*squeeze(U1(N,2,:))) - (s + 1i .*1.*omega1) .* M\L * squeeze(U1(N,2,:))), 'facecolor', 'interp'); shading interp;

%% sanity check
% bn = elements.boundaryIdx;
% harmonicIdx = 2;
% Ux = (Gx * squeeze(U1(N,harmonicIdx,:))).';  % N x 1
% Uy = (Gy * squeeze(U1(N,harmonicIdx,:))).';
% normal_x = elements.boundaryNormals(:,1);
% normal_y = elements.boundaryNormals(:,2);
% gradNormal =  (Ux(bn).' .* normal_x + Uy(bn).' .* normal_y);
% excitationCheck = zeros(size(elements.points,1),1);
% excitationCheck(bn) = gamma*squeeze(U1(N,harmonicIdx,bn)) + gradNormal;
% % check the domain residual, expect the biggest error on the boundary
% domresidual = squeeze(-Mass*kappasq(:,harmonicIdx,1)).*squeeze(U1(N,harmonicIdx,:)) + (S * squeeze(U1(N,harmonicIdx,:))) + MassB * squeeze(U1(N,harmonicIdx,:)) - Mass * F1(harmonicIdx,:).';
% u1boundary(elements.boundaryIdx) = gamma.*u1(0,boundaryPointsSource(:,1), boundaryPointsSource(:,2)) + dot(squeeze(u1grad(0,boundaryPointsSource(:,1),boundaryPointsSource(:,2))).',boundaryPointsSourceNormals.').';

%% for the inverse problem we construct u0_1, u0_2, u0_3 (in this case as solution of our approx. scheme)
% at x0 = (s^0, b^0, 0), where s^0, b^0 are space constant functions
% for finer triangular meshes (high accuracy) these can be precomputed and
% stored to speed up computation
s0 = 500.*ones(size(s));
b0 = 150.*ones(size(b));
eta0 = zeros(size(eta));
kappasq0 = constructKappaReparameterized(elements, s0, b0, [omega1 omega2 omega1], N); % compute all the complex wave numbers needed
excitationsReferenceState = excitations;
% for i=1:size(elements.points,1)
%     u1Sampled(i) = u1F(elements.points(i,1), elements.points(i,2));
% end
% the boundary excitation we take from our exemplary function
[~, U0_1, F1] = solveWesterveltMultiLevelBoundaryExcitation(elements, omega1, beta, gamma, squeeze(kappasq0(:,:,1)), squeeze(excitationsReferenceState(:,:,1)), eta0, b0, nIter, N, 10^(-12));
[~, U0_2, F2] = solveWesterveltMultiLevelBoundaryExcitation(elements, omega2, beta, gamma, squeeze(kappasq0(:,:,2)), squeeze(excitationsReferenceState(:,:,2)), eta0, b0, nIter, N, 10^(-12));
[cN, U0_3, F3] = solveWesterveltMultiLevelBoundaryExcitation(elements, omega1, beta, gamma, squeeze(kappasq0(:,:,3)), squeeze(excitationsReferenceState(:,:,3)), eta0, b0, nIter, N, 10^(-12));

% there reference state in time x space
u0_1 = calcSolution(timeMesh, squeeze(U0_1(N,:,:)), omega1);
u0_2 = calcSolution(timeMesh, squeeze(U0_2(N,:,:)), omega2);
u0_3 = calcSolution(timeMesh, squeeze(U0_3(N,:,:)), omega1);

% just a quick sanity check --> the solutions MUST fulfill |\Delta u^0| \geq c > 0 and |u^0| \geq c > 0 a.e. in our time space cylinder
if min(min(u0_1)) <= 0 || min(min(u0_2)) <= 0 || min(min(u0_3)) <= 0
    error('One of the reference states is not bounded from below.');
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
end

u1b = calcRobinBoundary(elements, u1sampled, gamma, Gx, Gy);
u2b = calcRobinBoundary(elements, u2sampled, gamma, Gx, Gy);
u3b = calcRobinBoundary(elements, u3sampled, gamma, Gx, Gy);

% if these values are big, or do not decrease with a finer triangular mesh
% something is odd
l2boundaryValErroru0_1 = sqrt(sum(sum(abs(u1b(:,elements.boundaryIdx) - u0_1_rb(:,elements.boundaryIdx)).^2,2),1));
linfboundaryValErroru0_1 = sqrt(max(max(abs(u1b(:,elements.boundaryIdx) - u0_1_rb(:,elements.boundaryIdx)))));

l2boundaryValErroru0_2 = sqrt(sum(sum(abs(u2b(:,elements.boundaryIdx) - u0_2_rb(:,elements.boundaryIdx)).^2,2),1));
linfboundaryValErroru0_2 = sqrt(max(max(abs(u2b(:,elements.boundaryIdx) - u0_2_rb(:,elements.boundaryIdx)))));

l2boundaryValErroru0_3 = sqrt(sum(sum(abs(u3b(:,elements.boundaryIdx) - u0_3_rb(:,elements.boundaryIdx)).^2,2),1));
linfboundaryValErroru0_3 = sqrt(max(max(abs(u3b(:,elements.boundaryIdx) - u0_3_rb(:,elements.boundaryIdx)))));


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
u0_3sampledrecon = calcSolution(timeMesh, 2*u0sampled, omega1);


% another sanity check
if norm(norm(abs(u0_3sampledrecon - u3sampled),2),2) > 10e-8
    error('Fourier coefficients of reference state 3 do not match.');
end


%% 
% u0 is just a reference state - for the frozen Newton method
% x0 holds all the reference states and the respective initial values (all
% the same for each of the reference states)
x0.refState_1.u0 = u0sampled;
x0.refState_1.laplaceu0 = laplaceu0;
x0.refState_1.kappa0 = squeeze(kappasq0(:,:,1));

x0.refState_2.u0 = u0sampled;
x0.refState_2.laplaceu0 = laplaceu0;
x0.refState_2.kappa0 = squeeze(kappasq0(:,:,2));

x0.refState_3.u0 = 2*u0sampled;
x0.refState_3.laplaceu0 = 2*laplaceu0;
x0.refState_3.kappa0 = squeeze(kappasq0(:,:,3));

x0.refState_1.s0 = mean(s.*0.7).*ones(size(s));
x0.refState_1.b0 = mean(b.*0.7).*ones(size(s));
x0.refState_1.eta0 = eta0;

x0.refState_2.s0 = mean(s.*0.7).*ones(size(s));
x0.refState_2.b0 = mean(b.*0.7).*ones(size(s));
x0.refState_2.eta0 = eta0;

x0.refState_3.s0 = mean(s.*0.7).*ones(size(s));
x0.refState_3.b0 = mean(b.*0.7).*ones(size(s));
x0.refState_3.eta0 = eta0;

% the differences
dx.ds = s - x0.refState_1.s0;
dx.db = b -  x0.refState_1.b0;
dx.deta = eta -  x0.refState_1.eta0;
dx.excitation = zeros(size(elements.points,1), N);


[~, DU_1, DF_1] = solveLinearisedWesterveltMultiLevelBoundaryExcitation(elements, omega1, beta, gamma, x0.refState_1, dx, nIter, N, false);

[~, DU_2, DF_2] = solveLinearisedWesterveltMultiLevelBoundaryExcitation(elements, omega2, beta, gamma, x0.refState_2, dx, nIter, N, false);

[~, DU_3, DF_3] = solveLinearisedWesterveltMultiLevelBoundaryExcitation(elements, omega1, beta, gamma, x0.refState_3, dx, nIter, N, false);

% solution using (xn), we can compute the residual for each harmonic
% we need to defined the residual on the whole of \overline{\Omega}, the
% FEM solver demands this
residual_1 = zeros(size(squeeze(U1(N,:,:))));
residual_3 = zeros(size(squeeze(U1(N,:,:))));
residual_2 = zeros(size(squeeze(U1(N,:,:))));

% start with x0
xn = x0;;
% F_1 (x_n^\delta)
% reconstruct kappa for x_n as it depends on s and b
% kappasq = constructKappaReparameterized(elements, s, b, [omega1 omega2 omega1], N); % compute all the complex wave numbers needed

[~, Un_1, ~] = solveWesterveltMultiLevelBoundaryExcitation(elements, omega1, beta, gamma,  xn.refState_1.kappa0, squeeze(excitationsReferenceState(:,:,1)), xn.eta0, xn.b0, nIter, N, 10^(-12));

% r_1:= h^\delta_1 - F_1 (x_n^\delta)
residual_1(:,measurementPointsIdx) = squeeze(U1(N,:,measurementPointsIdx) - Un_1(N,:,measurementPointsIdx));
% the adjoint state (linearised PDE) is only driven by the observation
% difference, all the conjugation is handled by the function itself
%
[~, Uadj_1, Fadj_1, db_a, ds_a, deta_a] = solveAdjointLinearisedWesterveltMultiLevelBoundaryExcitation(elements, omega1, beta, gamma, x0.refState_1, residual_1, nIter, N);

db_a_t_1 = calcSolution(timeMesh, db_a, omega1);
db_a_int_1 = sum(db_a_t_1 .* timeMesh(2) - timeMesh(1),1); % TODO: use trapezoid rule

db_s_t_1 = calcSolution(timeMesh, ds_a, omega1);
db_s_int_1 = sum(db_s_t_1 .* timeMesh(2) - timeMesh(1),1); % TODO: use trapezoid rule

db_eta_t_1 = calcSolution(timeMesh, deta_a, omega1);
db_eta_int_1 = sum(db_eta_t_1 .* timeMesh(2) - timeMesh(1),1); % TODO: use trapezoid rule1

%%
[~, Un_2, ~] = solveWesterveltMultiLevelBoundaryExcitation(elements, omega2, beta, gamma, xn.refState_2.kappa0, squeeze(excitationsReferenceState(:,:,2)), xn.eta0, xn.b0, nIter, N, 10^(-12));
residual_2(:,measurementPointsIdx) = squeeze(U2(N,:,measurementPointsIdx) - Un_2(N,:,measurementPointsIdx));

[~, Uadj_2, Fadj_2, db_a, ds_a, deta_a] = solveAdjointLinearisedWesterveltMultiLevelBoundaryExcitation(elements, omega2, beta, gamma, x0.refState_2, residual_2, nIter, N);

db_a_t_2 = calcSolution(timeMesh, db_a, omega2);
db_a_int_2 = sum(db_a_t_2 .* timeMesh(2) - timeMesh(1),1); % TODO: use trapezoid rule

db_s_t_2 = calcSolution(timeMesh, ds_a, omega2);
db_s_int_2 = sum(db_s_t_2 .* timeMesh(2) - timeMesh(1),1); % TODO: use trapezoid rule

db_eta_t_2 = calcSolution(timeMesh, deta_a, omega2);
db_eta_int_2 = sum(db_eta_t_2 .* timeMesh(2) - timeMesh(1),1); % TODO: use trapezoid rule1


%%
[~, Un_3, ~] = solveWesterveltMultiLevelBoundaryExcitation(elements, omega1, beta, gamma, x0.refState_3.kappa0, squeeze(excitationsReferenceState(:,:,3)), xn.eta0, xn.b0, nIter, N, 10^(-12));
residual_3(:,measurementPointsIdx) = squeeze(U3(N,:,measurementPointsIdx) - Un_3(N,:,measurementPointsIdx));

[~, Uadj_3, Fadj_3, db_a, ds_a, deta_a] = solveAdjointLinearisedWesterveltMultiLevelBoundaryExcitation(elements, omega1, beta, gamma, x0.refState_3, residual_3, nIter, N);

db_a_t_3 = calcSolution(timeMesh, db_a, omega1);
db_a_int_3 = sum(db_a_t_3 .* timeMesh(2) - timeMesh(1),1); % TODO: use trapezoid rule

db_s_t_3 = calcSolution(timeMesh, ds_a, omega1);
db_s_int_3 = sum(db_s_t_3 .* timeMesh(2) - timeMesh(1),1); % TODO: use trapezoid rule

db_eta_t_3 = calcSolution(timeMesh, deta_a, omega1);
db_eta_int_3 = sum(db_eta_t_3 .* timeMesh(2) - timeMesh(1),1); % TODO: use trapezoid rule1

%% the names of the actual parameters are not yet nice, but it makes things a bit easier (TODO)
xn = x0; % start at x0
alpha = 1; % alpha0
q = 1/2;

newtonIterations = 10;

for newtonIter = 1:newtonIterations

    % in each Newton step we have to do a CG

    % A = K*K + P*P + alpha

    % A(xn)
    y = applyA(xn, alpha, elements, timeMesh, x0, beta, gamma, [omega1 omega2 omega1], nIter, N);

    %prepare rhs = A(xn) + K^*(h - F(xn)) + \alpha_n(x0 - xn)
    % F(xn)
    % do not forget to update kappa
    xn.refState_1.kappa0 = constructKappaReparameterized(elements, xn.refState_1.s0, xn.refState_1.b0, omega1, N);
    xn.refState_2.kappa0 = constructKappaReparameterized(elements, xn.refState_2.s0, xn.refState_2.b0, omega2, N);
    xn.refState_3.kappa0 = constructKappaReparameterized(elements, xn.refState_3.s0, xn.refState_3.b0, omega1, N);

    [~, Un_1, ~] = solveWesterveltMultiLevelBoundaryExcitation(elements, omega1, beta, gamma, xn.refState_1.kappa0, squeeze(excitationsReferenceState(:,:,1)), xn.refState_1.eta0, xn.refState_1.b0, nIter, N, 10^(-12));
    [~, Un_2, ~] = solveWesterveltMultiLevelBoundaryExcitation(elements, omega2, beta, gamma, xn.refState_2.kappa0, squeeze(excitationsReferenceState(:,:,2)), xn.refState_2.eta0, xn.refState_2.b0, nIter, N, 10^(-12));
    [~, Un_3, ~] = solveWesterveltMultiLevelBoundaryExcitation(elements, omega1, beta, gamma, xn.refState_3.kappa0, squeeze(excitationsReferenceState(:,:,3)), xn.refState_3.eta0, xn.refState_3.b0, nIter, N, 10^(-12));
    residual_1 = zeros(size(squeeze(Un_1(N,:,:))));
    residual_3 = zeros(size(squeeze(Un_1(N,:,:))));
    residual_2 = zeros(size(squeeze(Un_1(N,:,:))));
    residual_1(:,measurementPointsIdx) = squeeze(U1(N,:,measurementPointsIdx) - Un_1(N,:,measurementPointsIdx));
    residual_2(:,measurementPointsIdx) = squeeze(U2(N,:,measurementPointsIdx) - Un_2(N,:,measurementPointsIdx));
    residual_3(:,measurementPointsIdx) = squeeze(U3(N,:,measurementPointsIdx) - Un_3(N,:,measurementPointsIdx));

    %K^*(h  - F(xn))
    [~, Uadj_1, Fadj_1, db_a, ds_a, deta_a] = solveAdjointLinearisedWesterveltMultiLevelBoundaryExcitation(elements, omega1, beta, gamma, x0.refState_1, residual_1, nIter, N);
    [db_1, ds_1, deta_1] = calcAdjointParameterStates(db_a, ds_a, deta_a, timeMesh, omega1);

    [~, Uadj_2, Fadj_2, db_a, ds_a, deta_a] = solveAdjointLinearisedWesterveltMultiLevelBoundaryExcitation(elements, omega2, beta, gamma, x0.refState_2, residual_2, nIter, N);
    [db_2, ds_2, deta_2] = calcAdjointParameterStates(db_a, ds_a, deta_a, timeMesh, omega2);

    [~, Uadj_3, Fadj_3, db_a, ds_a, deta_a] = solveAdjointLinearisedWesterveltMultiLevelBoundaryExcitation(elements, omega1, beta, gamma, x0.refState_3, residual_3, nIter, N);
    [db_3, ds_3, deta_3] = calcAdjointParameterStates(db_a, ds_a, deta_a, timeMesh, omega1);
    rhs = xn;
    rhs.refState_1.eta0   = y.refState_1.eta0   + deta_1 + alpha.* (x0.refState_1.eta0 - xn.refState_1.eta0);
    rhs.refState_1.s0     = y.refState_1.s0     + ds_1   + alpha.* (x0.refState_1.s0 - xn.refState_1.s0);
    rhs.refState_1.b0     = y.refState_1.b0     + db_1   + alpha.* (x0.refState_1.b0 - xn.refState_1.b0);

    rhs.refState_2.eta0   = y.refState_2.eta0   + deta_2 + alpha.* (x0.refState_2.eta0 - xn.refState_2.eta0);
    rhs.refState_2.s0     = y.refState_2.s0     + ds_2   + alpha.* (x0.refState_2.s0 - xn.refState_2.s0);
    rhs.refState_2.b0     = y.refState_2.b0     + db_2   + alpha.* (x0.refState_2.b0 - xn.refState_2.b0);

    rhs.refState_3.eta0   = y.refState_3.eta0  + deta_3  + alpha.* (x0.refState_3.eta0 - xn.refState_3.eta0);
    rhs.refState_3.s0     = y.refState_3.s0    + ds_3    + alpha.* (x0.refState_3.s0 - xn.refState_3.s0);
    rhs.refState_3.b0     = y.refState_3.b0    + db_3    + alpha.* (x0.refState_3.b0 - xn.refState_3.b0);

    % now we need to solve Az = rhs
    CGIterations = 40;
    z = xn;
    res = minusParameters(rhs, applyA(z, alpha, elements, timeMesh, x0, beta, gamma, [omega1 omega2 omega1], nIter, N)); %
    pk = res;
    stopres = zeros(CGIterations,1);
    for iter = 1:CGIterations
        Apk = applyA(pk, alpha, elements, timeMesh, x0, beta, gamma, [omega1 omega2 omega1], nIter, N);
        rr = calcInnerProductParameters(res,res, elements);
        d(iter) = rr / (calcInnerProductParameters(pk,Apk ,elements));
        z = addParameters(z, scalarMulParameters(d(iter), pk));
        resNew = minusParameters(res, scalarMulParameters(d(iter), Apk));
        rrN = calcInnerProductParameters(resNew, resNew, elements);
        stopres(iter) =rrN;

        if stopres(iter) < 10e-20
            break;
        end

        betak(iter) = rrN/rr;

        pk = addParameters(resNew, scalarMulParameters(betak(iter), pk));
        res = resNew;
    end

    xn = z;
    alpha = alpha*q;
end

%%
point = [0.0;0.05];

[v,idx] = min(sum((elements.points - point(:)').^2,2)); 

node = [elements.points(idx,1);elements.points(idx,2)];

T = 1/f2;
% sampling frequency in time
Fs = 1/T * 2 * (N);
omega = omega2;
U = squeeze(U0_2(N,:,:));

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
