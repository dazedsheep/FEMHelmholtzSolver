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

du = calcSolution(timeMesh.timeMesh1, squeeze(DU(N,:,:)), omega1);

err = u1s + du - u1s;
if max(max(abs(err))) > 10e-15
    error("DF(x_0)(0) is not zero");
end

% now with a small pertubation
perturbationCoeff = 0.99;
perturbed_s = s.*perturbationCoeff;
perturbed_b = b.*perturbationCoeff;
perturbed_eta = eta.*perturbationCoeff;

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
du = calcSolution(timeMesh.timeMesh1, squeeze(DU(N,:,:)), omega1);
[~, ldx.ds] = integrate_fun_trimesh(elements.opoints, elements.otri, abs(dx.ds).^2.');
[~, ldx.db] = integrate_fun_trimesh(elements.opoints, elements.otri, abs(dx.db).^2.');
[~, ldx.deta] = integrate_fun_trimesh(elements.opoints, elements.otri, abs(dx.deta).^2.');

nm = sqrt(ldx.ds + ldx.db + ldx.deta);
% for very small pertubations we should have F(x) \approx F(x0)
diffU = (UPert(N,:,:) - DU(N,:,:) - U1(N,:,:));
diffU_time = calcSolution(timeMesh.timeMesh1, squeeze(diffU), omega1);
%figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2),trapz(timeMesh.timeMesh1, abs(calcSolution(timeMesh.timeMesh1,squeeze(UPert(N,:,:) - U1(N,:,:)), omega1)).^2,1), 'facecolor', 'interp'); shading interp;
%figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2),trapz(timeMesh.timeMesh1, abs(du).^2,1), 'facecolor', 'interp'); shading interp;
%du_diff = calcSolution(timeMesh.timeMesh1,squeeze(U1(N,:,:) - UPert(N,:,:)), omega1);
%figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2), du_diff(1,:), 'facecolor', 'interp'); shading interp;
%figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2), du(1,:), 'facecolor', 'interp'); shading interp;

for i = 1:size(timeMesh.timeMesh1,2)
    [~, a(i)] = integrate_fun_trimesh(elements.opoints, elements.otri, abs(diffU_time(1,:)).^2);
end
dn = trapz(timeMesh.timeMesh1,a);

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
du = calcSolution(timeMesh.timeMesh1, squeeze(DU(N,:,:)), omega1);
[~, ldx.ds] = integrate_fun_trimesh(elements.opoints, elements.otri, abs(dx.ds).^2.');
[~, ldx.db] = integrate_fun_trimesh(elements.opoints, elements.otri, abs(dx.db).^2.');
[~, ldx.deta] = integrate_fun_trimesh(elements.opoints, elements.otri, abs(dx.deta).^2.');

nm = sqrt(ldx.ds + ldx.db + ldx.deta);
% for very small pertubations we should have F(x) \approx F(x0)
diffU = (UPert(N,:,:) - DU(N,:,:) - U1(N,:,:));
diffU_time = calcSolution(timeMesh.timeMesh1, squeeze(diffU), omega1);
for i = 1:size(timeMesh.timeMesh1,2)
    [~, a(i)] = integrate_fun_trimesh(elements.opoints, elements.otri, abs(diffU_time(1,:)).^2);
end
dn = trapz(timeMesh.timeMesh1,a);
harmonics = 0:(N-1);
lapu0 = squeeze(kappaPerturbed(:,:,1)).'.*squeeze(UPert(N,:,:));
calcLap0sol = calcSolution(timeMesh.timeMesh1, -harmonics.'.^2 .* lapu0, omega1);
% calc also the adjoint state
residual_1 = zeros(N,size(elements.points,1));
residual_1(:,elements.measurementPointsIdx) = squeeze(DU(N,:,elements.measurementPointsIdx));
[~, Uadj_1, ~] = solveAdjointLinearisedWesterveltMultiLevelBoundaryExcitation(elements, omega1, beta, gamma, xs.refState_1, residual_1, nIter, N);
uadj = calcSolution(timeMesh.timeMesh1, squeeze(Uadj_1(N,:,:)), omega1);
ads = trapz(timeMesh.timeMesh1, uadj.*(calcLap0sol));
[~,int_ds] = integrate_fun_trimesh(elements.opoints, elements.otri, dx.ds.*ads);


%% the names of the actual parameters are not yet nice, but it makes things a bit easier (TODO)
% u0 is just a reference state - for the frozen Newton method
% x0 holds all the reference states and the respective initial values (all
% the same for each of the reference states)
kappasq0 = constructKappaReparameterized(elements, s0, b0, [omega1 omega2 omega1], N); % compute all the complex wave numbers needed
excitationsReferenceState = excitations;
alpha = 1;

x0.refState_1.u0 = u0sampled;
x0.refState_1.laplaceu0 = laplaceu0;
x0.refState_1.kappa0 = squeeze(kappasq0(:,:,1));

x0.refState_2.u0 = u0sampled;
x0.refState_2.laplaceu0 = laplaceu0;
x0.refState_2.kappa0 = squeeze(kappasq0(:,:,2));

x0.refState_3.u0 = u3Amplitude*u0sampled;
x0.refState_3.laplaceu0 = u3Amplitude*laplaceu0;
x0.refState_3.kappa0 = squeeze(kappasq0(:,:,3));

x0.refState_1.s0 = s0;
x0.refState_1.b0 = b0;
x0.refState_1.eta0 = eta0;

x0.refState_2.s0 = s0;
x0.refState_2.b0 = b0;
x0.refState_2.eta0 = eta0;

x0.refState_3.s0 = s0;
x0.refState_3.b0 = b0;
x0.refState_3.eta0 = eta0;

referenceStates.u1LaplacianSampled = u1LaplacianSampled;
referenceStates.u1ttSampled = u1ttSampled;
referenceStates.u1sqttSampled = u1sqttSampled;

referenceStates.u2LaplacianSampled = u2LaplacianSampled;
referenceStates.u2ttSampled = u2ttSampled;
referenceStates.u2sqttSampled = u2sqttSampled;

referenceStates.u3LaplacianSampled = u3LaplacianSampled;
referenceStates.u3ttSampled = u3ttSampled;
referenceStates.u3sqttSampled = u3sqttSampled;

%% we do some operator tests, to make sure the most principle things are fine
a = 100;
% the differences
dx1.ds = (s).*rand(1);
dx1.db = (b).*rand(1);
dx1.deta = (eta ).*rand(1);
dx1.excitation = zeros(size(elements.points,1), N);

dxa.ds = a.*dx1.ds;
dxa.db = a.*dx1.db;
dxa.deta = a.*dx1.deta;
dxa.excitation = zeros(size(elements.points,1), N);

dx2.ds = (s).*rand(1);
dx2.db = (b).*rand(1);
dx2.deta = (eta).*rand(1);
dx2.excitation = zeros(size(elements.points,1), N);

dxs.ds = dx1.ds + dx2.ds;
dxs.db = dx1.db + dx2.db;
dxs.deta = dx1.deta + dx2.deta;
dxs.excitation = zeros(size(elements.points,1), N);

[~, du1, ~] = solveLinearisedWesterveltMultiLevelBoundaryExcitation(elements, omega1, beta, gamma, x0.refState_1, dx1, nIter, N, false);
[~, du2, ~] = solveLinearisedWesterveltMultiLevelBoundaryExcitation(elements, omega1, beta, gamma, x0.refState_1, dx2, nIter, N, false);
[~, dus, ~] = solveLinearisedWesterveltMultiLevelBoundaryExcitation(elements, omega1, beta, gamma, x0.refState_1, dxs, nIter, N, false);
[~, dua, ~] = solveLinearisedWesterveltMultiLevelBoundaryExcitation(elements, omega1, beta, gamma, x0.refState_1, dxa, nIter, N, false);

% linear
diff1 = norm(norm(squeeze(du1(N,:,:) + du2(N,:,:) - dus(N,:,:))));
if diff1 > 10e-12
    error("Frechet derivative is not linear!");
end

% homogenity
diff2 = norm(norm(squeeze(dua(N,:,:) - a.*du1(N,:,:))));
if diff2 > 10e-10
    error("Frechet derivative is not linear (homogenity failed)!");
end


%% linearity test for operator A
a = 100.*rand(1);
x1 = x0;
x2 = x0;

x1.refState_1.s0 = s0.*rand(1);
x1.refState_1.b0 = b0.*rand(1);
x1.refState_1.eta0 = eta.*rand(1);

x1.refState_2.s0 = s0.*rand(1);
x1.refState_2.b0 = b0.*rand(1);
x1.refState_2.eta0 = eta.*rand(1);

x1.refState_3.s0 = s0.*rand(1);
x1.refState_3.b0 = b0.*rand(1);
x1.refState_3.eta0 = eta.*rand(1);

x2.refState_1.s0 = s0.*rand(1);
x2.refState_1.b0 = b0.*rand(1);
x2.refState_1.eta0 = eta.*rand(1);

x2.refState_2.s0 = s0.*rand(1);
x2.refState_2.b0 = b0.*rand(1);
x2.refState_2.eta0 = eta.*rand(1);

x2.refState_3.s0 = s0.*rand(1);
x2.refState_3.b0 = b0.*rand(1);
x2.refState_3.eta0 = eta.*rand(1);

Ax1 = applyA(x1, alpha, elements, timeMesh, x0, beta, gamma, [omega1 omega2 omega1], nIter, N, referenceStates, useSolutionAsLinPoint);
Ax2 = applyA(x2, alpha, elements, timeMesh, x0, beta, gamma, [omega1 omega2 omega1], nIter, N, referenceStates, useSolutionAsLinPoint);
Ax3 = applyA(addParameters(x1,x2), alpha, elements, timeMesh, x0, beta, gamma, [omega1 omega2 omega1], nIter, N, referenceStates, useSolutionAsLinPoint);
Axh = applyA(scalarMulParameters(a,x2), alpha, elements, timeMesh, x0, beta, gamma, [omega1 omega2 omega1], nIter, N, referenceStates, useSolutionAsLinPoint);

diff1 = minusParameters(addParameters(Ax1, Ax2), Ax3);

if calcInnerProductParameters(diff1,diff1,elements) > 10e-12
    error("Operator A is not linear!");
end

diff2 = minusParameters(scalarMulParameters(a,Ax2), Axh);

if calcInnerProductParameters(diff2,diff2,elements) > 10e-12
    error("Operator A is not linear (homogenity failed)!");
end

%% symmetric and psd test for a part of the operator A
a = 2;
dx1.ds = s;
dx1.db = b;
dx1.deta = eta;
dx1.excitation = zeros(size(elements.points,1), N);
residual_1 = zeros(N,size(elements.points,1));

[~, DU_1, ~] = solveLinearisedWesterveltMultiLevelBoundaryExcitation(elements, omega1, beta, gamma, x0.refState_1, dx1, nIter, N, useSolutionAsLinPoint);
residual_1(:,elements.measurementPointsIdx) = squeeze(DU_1(N,:,elements.measurementPointsIdx));

[~, Uadj_1, ~] = solveAdjointLinearisedWesterveltMultiLevelBoundaryExcitation(elements, omega1, beta, gamma, x0.refState_1, residual_1, nIter, N);
[db_1, ds_1, deta_1] = calcAdjointStates((squeeze(Uadj_1(N,:,:))), omega1, timeMesh, referenceStates.u1LaplacianSampled, referenceStates.u1ttSampled, referenceStates.u1sqttSampled);

dx1.ds = a.*s;
dx1.db = b;
dx1.deta = a.*eta;


[~, DU_2, ~] = solveLinearisedWesterveltMultiLevelBoundaryExcitation(elements, omega1, beta, gamma, x0.refState_1, dx1, nIter, N, useSolutionAsLinPoint);
residual_1(:,elements.measurementPointsIdx) = squeeze(DU_2(N,:,elements.measurementPointsIdx));
[~, Uadj_2, ~] = solveAdjointLinearisedWesterveltMultiLevelBoundaryExcitation(elements, omega1, beta, gamma, x0.refState_1, residual_1, nIter, N);
[db_2, ds_2, deta_2] = calcAdjointStates((squeeze(Uadj_2(N,:,:))), omega1, timeMesh, referenceStates.u1LaplacianSampled, referenceStates.u1ttSampled, referenceStates.u1sqttSampled);


%% symmetric test for the operator A (a very simple one)

x1 = x0;
x2 = x0;
a = 10;
alpha = 1;
x1.refState_1.s0 = s.*rand(1).*a;
x1.refState_1.b0 = b.*rand(1).*a;
x1.refState_1.eta0 = eta.*rand(1).*a;

x1.refState_2.s0 = s.*rand(1).*a;
x1.refState_2.b0 = b.*rand(1).*a;
x1.refState_2.eta0 = eta.*rand(1).*a;

x1.refState_3.s0 = s.*rand(1).*a;
x1.refState_3.b0 = b.*rand(1).*a;
x1.refState_3.eta0 = eta.*rand(1).*a;

x2.refState_1.s0 = s.*rand(1).*a;
x2.refState_1.b0 = b.*rand(1).*a;
x2.refState_1.eta0 = eta.*rand(1).*a;

x2.refState_2.s0 = s.*rand(1).*a;
x2.refState_2.b0 = b.*rand(1).*a;
x2.refState_2.eta0 = eta.*rand(1).*a;

x2.refState_3.s0 = s.*rand(1).*a;
x2.refState_3.b0 = b.*rand(1).*a;
x2.refState_3.eta0 = eta.*rand(1).*a;

Ax1 = applyA(x1, alpha, elements, timeMesh, x0, beta, gamma, [omega1 omega2 omega3], nIter, N, referenceStates, useSolutionAsLinPoint);
Ax2 = applyA(x2, alpha, elements, timeMesh, x0, beta, gamma, [omega1 omega2 omega3], nIter, N, referenceStates, useSolutionAsLinPoint);

% <Ax1, x2> = <x1, A*x2> = <x1, Ax2>
innerP1 = calcInnerProductParameters(Ax1, x2 ,elements);
innerP2 = calcInnerProductParameters(x1, Ax2 ,elements);

if innerP1 < 0
    error("A is not p.s.d!");
end

diff = abs(innerP1 - innerP2);
if diff > 10e-3
    warning("Operator A is not self-adjoint");
end

%% do some sanity checks for the adjoint
dx.ds = s - x0.refState_1.s0;
dx.db = b - b;
dx.deta = eta  - x0.refState_1.eta0;
dx.excitation = zeros(size(elements.points,1), N);

[~, DU, ~] = solveLinearisedWesterveltMultiLevelBoundaryExcitation(elements, omega1, beta, gamma, x0.refState_1, dx, nIter, N, true);
du = calcSolution(timeMesh.timeMesh1, squeeze(DU(N,:,:)), omega1);

residual_1(:,elements.measurementPointsIdx) = squeeze(DU(N,:,elements.measurementPointsIdx));
[~, Uadj_1, ~] = solveAdjointLinearisedWesterveltMultiLevelBoundaryExcitation(elements, omega1, beta, gamma, x0.refState_1, residual_1, nIter, N);
uadj = calcSolution(timeMesh.timeMesh1, squeeze((Uadj_1(N,:,:))), omega1);
ads = trapz(timeMesh.timeMesh1, uadj.*referenceStates.u1LaplacianSampled);
adb = trapz(timeMesh.timeMesh1, uadj.*referenceStates.u1ttSampled);
adeta = trapz(timeMesh.timeMesh1, uadj.*referenceStates.u1sqttSampled);
[adb1,ads1, adeta1] = calcAdjointStates(squeeze((Uadj_1(N,:,:))), omega1, timeMesh.timeMesh1, referenceStates.u1LaplacianSampled, referenceStates.u1ttSampled, referenceStates.u1sqttSampled);
[~,int_ds] = integrate_fun_trimesh(elements.opoints, elements.otri, dx.ds.*ads);
[~,int_db] = integrate_fun_trimesh(elements.opoints, elements.otri, dx.db.*adb);
[~,int_deta] = integrate_fun_trimesh(elements.opoints, elements.otri, dx.deta.*adeta);
zz = -int_ds + int_db - int_deta;
