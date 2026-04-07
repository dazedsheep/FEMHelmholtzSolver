function [i, u, F] = solveWesterveltMultiLevelBoundaryExcitation(elements, omega, beta, gamma, kappa, excitation, f, b, nIterations, nHarmonics, threshold)

n = size(elements.points,1);
N = nIterations;
u = zeros(N,N,n);
F = zeros(N, n);

for i=1:N
    for j=0:min((i-1),nHarmonics-1)
        p_m = zeros(1,n);
        % the first iteration has just the excitation on the right hand side
        % index 1 is the zero-th solution
        if i>1 && j>0
            % ---------- First sum ----------
            % sum_{l=0}^j u_l * u_{j-l}
            for l = 0:j
                p_m = p_m + ...
                    squeeze(u(i-1,l+1,:)).' .* ...
                    squeeze(u(i-1,(j-l)+1,:)).';
            end

            % ---------- Second sum ----------
            % 2 * sum_{r=0}^{N-1-j} conj(u_r) * u_{r+j}
             for r = j:2:(2*(nHarmonics-1) - j)
                minusidx = (r-j)/2;
                plusidx = (r+j)/2;
                p_m = p_m + 2 * ...
                    conj(squeeze(u(i-1,minusidx+1,:)).') .* ...
                    squeeze(u(i-1,plusidx+1,:)).';
            end

        end
        
        F(j+1,:) = -j^2.*kappa(:,j+1).*1./(2.*b).*f.*p_m.' + excitation(:, j+1);
        F(j+1,elements.boundaryIdx) = 0;

        u(i,j+1,:) = solveHelmholtzCondensedC(elements, j*omega, gamma, j^2.*kappa(:,j+1), beta, F(j+1,:).', excitation(:,j+1), n,  elements.K, elements.rowK, elements.colK, elements.M_t, elements.tBM);

    end

end


end