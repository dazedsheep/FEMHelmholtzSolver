function [xn] = frozenNewtonMethod(elements, timeMesh, U0_1, U0_F1, U0_2, U0_F2, U0_3, U0_F3, s0, b0, eta0, beta, gamma, measU1, measU2, measU3, omega1, omega2, omega3, excitations, useSolutionAsLinPoint, nIterations, nHarmonics, newtonIterations, NewtonTol, CGTol)

N = nHarmonics;
nIter = nIterations;

excitationsReferenceState = excitations;

alpha = 1; % alpha0
q = 1/2;

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

    x0.refState_3.u0 = u3Amplitude*u0sampled;
    x0.refState_3.laplaceu0 = u3Amplitude*laplaceu0;
    x0.refState_3.kappa0 = squeeze(kappasq0(:,:,3));

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

CGIterations = 50;

xn = x0; % start at x0
residue = ones(newtonIterations,4);

for newtonIter = 1:newtonIterations

    % in each Newton step we have to do a CG

    % A = K*K + P*P + alpha
    %xn = alignParameters(xn);
    % A(xn)
    y = applyA(xn, alpha, elements, timeMesh, x0, beta, gamma, [omega1 omega2 omega3], nIter, N, referenceStates, useSolutionAsLinPoint);

    %prepare rhs = A(xn) + K^*(h - F(xn)) + \alpha_n(x0 - xn)
    % F(xn)
    % do not forget to update kappa

    kappa1 = constructKappaReparameterized(elements, xn.refState_1.s0, xn.refState_1.b0, omega1, N);
    kappa2 = constructKappaReparameterized(elements, xn.refState_2.s0, xn.refState_2.b0, omega2, N);
    kappa3 = constructKappaReparameterized(elements, xn.refState_3.s0, xn.refState_3.b0, omega3, N);


    [~, Un_1, ~] = solveWesterveltMultiLevelBoundaryExcitation(elements, omega1, beta, gamma, kappa1, squeeze(excitationsReferenceState(:,:,1)), xn.refState_1.eta0, xn.refState_1.b0, nIter, N, 10^(-12));
    [~, Un_2, ~] = solveWesterveltMultiLevelBoundaryExcitation(elements, omega2, beta, gamma, kappa2, squeeze(excitationsReferenceState(:,:,2)), xn.refState_2.eta0, xn.refState_2.b0, nIter, N, 10^(-12));
    [~, Un_3, ~] = solveWesterveltMultiLevelBoundaryExcitation(elements, omega3, beta, gamma, kappa3, squeeze(excitationsReferenceState(:,:,3)), xn.refState_3.eta0, xn.refState_3.b0, nIter, N, 10^(-12));
    
    residual_1 = zeros(size(squeeze(Un_1(N,:,:))));
    residual_2 = zeros(size(squeeze(Un_2(N,:,:))));
    residual_3 = zeros(size(squeeze(Un_3(N,:,:))));
    
    residual_1(:,elements.measurementPointsIdx) = measU1(:,elements.measurementPointsIdx) - squeeze(Un_1(N,:,elements.measurementPointsIdx));
    residual_2(:,elements.measurementPointsIdx) = measU2(:,elements.measurementPointsIdx) - squeeze(Un_2(N,:,elements.measurementPointsIdx));
    residual_3(:,elements.measurementPointsIdx) = measU3(:,elements.measurementPointsIdx) - squeeze(Un_3(N,:,elements.measurementPointsIdx));

    % residual in L^2(\Sigma)
    [~, residue(newtonIter, 1)] = integrate_fun_trimesh(elements.opoints, elements.otri, sum(abs(residual_1).^2,1));
    [~, residue(newtonIter, 2)] = integrate_fun_trimesh(elements.opoints, elements.otri, sum(abs(residual_2).^2,1));
    [~, residue(newtonIter, 3)] = integrate_fun_trimesh(elements.opoints, elements.otri, sum(abs(residual_3).^2,1));
    residue(newtonIter, 4) = sqrt(residue(newtonIter, 1) + residue(newtonIter, 2) + residue(newtonIter, 3));
    fprintf('Current residue: %e\n',residue(newtonIter, 4));
    if residue(newtonIter, 4) < NewtonTol
        break;
    end

    %K^*(h  - F(xn))
    [~, Uadj_1, Fadj_1] = solveAdjointLinearisedWesterveltMultiLevelBoundaryExcitation(elements, omega1, beta, gamma, x0.refState_1, residual_1, nIter, N);
    [db_1, ds_1, deta_1] = calcAdjointStates((squeeze(Uadj_1(N,:,:))), omega1, timeMesh, referenceStates.u1LaplacianSampled, referenceStates.u1ttSampled, referenceStates.u1sqttSampled);

    [~, Uadj_2, Fadj_2] = solveAdjointLinearisedWesterveltMultiLevelBoundaryExcitation(elements, omega2, beta, gamma, x0.refState_2, residual_2, nIter, N);
    [db_2, ds_2, deta_2] = calcAdjointStates((squeeze(Uadj_2(N,:,:))), omega2, timeMesh, referenceStates.u2LaplacianSampled, referenceStates.u2ttSampled, referenceStates.u2sqttSampled);

    [~, Uadj_3, Fadj_3] = solveAdjointLinearisedWesterveltMultiLevelBoundaryExcitation(elements, omega3, beta, gamma, x0.refState_3, residual_3, nIter, N);
    [db_3, ds_3, deta_3] = calcAdjointStates((squeeze(Uadj_3(N,:,:))), omega3, timeMesh, referenceStates.u3LaplacianSampled, referenceStates.u3ttSampled, referenceStates.u3sqttSampled);

    rhs = xn;

    rhs.refState_1.eta0   = y.refState_1.eta0   + deta_1 + alpha.* (x0.refState_1.eta0  - xn.refState_1.eta0);
    rhs.refState_1.s0     = y.refState_1.s0     + ds_1   + alpha.* (x0.refState_1.s0    - xn.refState_1.s0);
    rhs.refState_1.b0     = y.refState_1.b0     + db_1   + alpha.* (x0.refState_1.b0    - xn.refState_1.b0);

    rhs.refState_2.eta0   = y.refState_2.eta0   + deta_2 + alpha.* (x0.refState_2.eta0  - xn.refState_2.eta0);
    rhs.refState_2.s0     = y.refState_2.s0     + ds_2   + alpha.* (x0.refState_2.s0    - xn.refState_2.s0);
    rhs.refState_2.b0     = y.refState_2.b0     + db_2   + alpha.* (x0.refState_2.b0    - xn.refState_2.b0);

    rhs.refState_3.eta0   = y.refState_3.eta0  + deta_3  + alpha.* (x0.refState_3.eta0  - xn.refState_3.eta0);
    rhs.refState_3.s0     = y.refState_3.s0    + ds_3    + alpha.* (x0.refState_3.s0    - xn.refState_3.s0);
    rhs.refState_3.b0     = y.refState_3.b0    + db_3    + alpha.* (x0.refState_3.b0    - xn.refState_3.b0);

    % now we need to solve Az = rhs
    z = xn;
    res = minusParameters(rhs, applyA(z, alpha, elements, timeMesh, x0, beta, gamma, [omega1 omega2 omega3], nIter, N, referenceStates, useSolutionAsLinPoint));
    pk = res;
    stopres = zeros(CGIterations,1);
    betak = zeros(CGIterations,1);
    d = betak;
    for iter = 1:CGIterations
        Apk = applyA(pk, alpha, elements, timeMesh, x0, beta, gamma, [omega1 omega2 omega3], nIter, N, referenceStates, useSolutionAsLinPoint);
        rr = calcInnerProductParameters(res,res, elements);
        d(iter) = rr / (calcInnerProductParameters(pk, Apk ,elements));
        z = addParameters(z, scalarMulParameters(d(iter), pk));
        % the following line may propagate numerical errors
        resNew = minusParameters(res, scalarMulParameters(d(iter), Apk));
        rrN = calcInnerProductParameters(resNew, resNew, elements);
        stopres(iter) = rrN;

        if stopres(iter) < CGTol
            break;
        end

        betak(iter) = rrN/rr;

        pk = addParameters(resNew, scalarMulParameters(betak(iter), pk));
        res = resNew;
    end

    xn = z;
    alpha = alpha*q;
end


end

