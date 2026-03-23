function [y] = applyA(xn, alpha, elements, timeMesh, x0, beta, gamma, omega, nIterations, nHarmonics, referenceStates, useSolutionAsLinPoint)
% this is the implementation of the operator A*A + P*P + alpha_n
N=nHarmonics;
nIter = nIterations;

residual_1 = zeros(N,size(elements.points,1));
residual_2 = zeros(N,size(elements.points,1));
residual_3 = zeros(N,size(elements.points,1));

% j = 1, K(x_n)
dx.deta = xn.refState_1.eta0;
dx.ds   = xn.refState_1.s0;
dx.db   = xn.refState_1.b0;
dx.excitation = zeros(size(elements.points,1), N);

[~, DU_1, ~] = solveLinearisedWesterveltMultiLevelBoundaryExcitation(elements, omega(1), beta, gamma, x0.refState_1, dx, nIter, N, useSolutionAsLinPoint);
residual_1(:,elements.measurementPointsIdx) = squeeze(DU_1(N,:,elements.measurementPointsIdx));
% the adjoint state (linearised PDE) is only driven by the observation
% difference, all the conjugation is handled by the function itself
%
% j = 1, K^*(K(x_n))
[~, Uadj_1, ~] = solveAdjointLinearisedWesterveltMultiLevelBoundaryExcitation(elements, omega(1), beta, gamma, x0.refState_1, residual_1, nIter, N);

[db_int_1, ds_int_1, deta_int_1] = calcAdjointStates((squeeze(Uadj_1(N,:,:))), omega(1), timeMesh, referenceStates.u1LaplacianSampled, referenceStates.u1ttSampled, referenceStates.u1sqttSampled);

%%
% j = 2, K(x_n)
dx.deta = xn.refState_2.eta0;
dx.ds   = xn.refState_2.s0;
dx.db   = xn.refState_2.b0;
[~, DU_2, ~] = solveLinearisedWesterveltMultiLevelBoundaryExcitation(elements, omega(2), beta, gamma, x0.refState_2, dx, nIter, N, useSolutionAsLinPoint);
residual_2(:,elements.measurementPointsIdx) = squeeze(DU_2(N,:,elements.measurementPointsIdx));
% j = 2, K^*(K(x_n))
[~, Uadj_2, ~] = solveAdjointLinearisedWesterveltMultiLevelBoundaryExcitation(elements, omega(2), beta, gamma, x0.refState_2, residual_2, nIter, N);

[db_int_2, ds_int_2, deta_int_2] = calcAdjointStates((squeeze(Uadj_2(N,:,:))), omega(2), timeMesh, referenceStates.u2LaplacianSampled, referenceStates.u2ttSampled, referenceStates.u2sqttSampled);

%%
% j = 3, K(x_n)
dx.deta = xn.refState_3.eta0;
dx.ds = xn.refState_3.s0;
dx.db = xn.refState_3.b0;
[~, DU_3, ~] = solveLinearisedWesterveltMultiLevelBoundaryExcitation(elements, omega(3), beta, gamma, x0.refState_3, dx, nIter, N, useSolutionAsLinPoint);
residual_3(:,elements.measurementPointsIdx) = squeeze(DU_3(N,:,elements.measurementPointsIdx));
% j = 3, K^*(K(x_n))
[~, Uadj_3, ~] = solveAdjointLinearisedWesterveltMultiLevelBoundaryExcitation(elements, omega(3), beta, gamma, x0.refState_3, residual_3, nIter, N);

[db_int_3, ds_int_3, deta_int_3] = calcAdjointStates((squeeze(Uadj_3(N,:,:))), omega(3), timeMesh, referenceStates.u3LaplacianSampled, referenceStates.u3ttSampled, referenceStates.u3sqttSampled);

% we do not need P*P
%%
% K^*(K(x)) + \alpha*xn
y = xn;
y.refState_1.eta0   = deta_int_1 + alpha.* xn.refState_1.eta0;
y.refState_1.s0     = ds_int_1   + alpha.* xn.refState_1.s0;
y.refState_1.b0     = db_int_1   + alpha.* xn.refState_1.b0;

y.refState_2.eta0   = deta_int_2 + alpha.* xn.refState_2.eta0;
y.refState_2.s0     = ds_int_2   + alpha.* xn.refState_2.s0;
y.refState_2.b0     = db_int_2   + alpha.* xn.refState_2.b0;

y.refState_3.eta0   = deta_int_3 + alpha.* xn.refState_3.eta0;
y.refState_3.s0     = ds_int_3   + alpha.* xn.refState_3.s0;
y.refState_3.b0     = db_int_3   + alpha.* xn.refState_3.b0;


end

