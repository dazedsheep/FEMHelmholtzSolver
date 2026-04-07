function [i, u, F] = solveAdjointLinearisedWesterveltMultiLevelBoundaryExcitation(elements, omega, beta, gamma, x0, yobs, nIterations, nHarmonics)
% yobs has to be defined on the hole of \overline{\Omega}, just fill stuff
% with 0
n = size(elements.points,1);
N = nIterations;
h = zeros(n,1);
u = zeros(N,N,n);
F = zeros(N, n);
for i=1:N
    for j=0:min((i-1),nHarmonics-1)
        p_m_eta = zeros(1,n);
        
        % the first iteration has just the excitation on the right hand side
        % index 1 is the zero-th solution
        if i>1 && j>0
            % for the linearised operator the right hand side is a bit more
            % complex

            % ---------- First sum ----------
            % sum_{l=1}^j u_l * p_{j-l} * j^2 * \omega^2
            for l = 0:j
                p_m_eta = p_m_eta + l.^2.*omega^2.* conj( ...
                    squeeze(x0.u0(l+1,:))) .* ...
                    squeeze(u(i-1,(j-l)+1,:)).';
            end

            % ---------- Second sum ----------
             for r = j:2:(2*(nHarmonics-1) - j)
                minusidx = (r-j)/2;
                plusidx = (r+j)/2;
                p_m_eta = p_m_eta + r.^2.* (...
                    x0.u0(minusidx+1,:).* ...
                    conj(squeeze(squeeze(u(i-1,(plusidx)+1,:)).')).*(plusidx).^2 + ...
                    (squeeze(u(i-1,minusidx+1,:)).').* ...
                    conj(squeeze(x0.u0(plusidx+1,:))).*minusidx.^2);
            end

        end
        cs = 1./(x0.s0.' - 1i.*j.*omega);
        F(j+1,:) = -x0.eta0.'.*conj(p_m_eta).*cs; % et0 part (the only one in the adjoint)
        F(j+1,elements.boundaryIdx) = 0;

        u(i,j+1,:) = solveHelmholtzCondensedC(elements, j*omega, gamma, conj(j^2.*x0.kappa0(:,j+1)), beta, F(j+1,:).', cs.'.*yobs(j+1,:).', n, elements.K, elements.rowK, elements.colK, elements.M_t, elements.tBM);
    end

end


end