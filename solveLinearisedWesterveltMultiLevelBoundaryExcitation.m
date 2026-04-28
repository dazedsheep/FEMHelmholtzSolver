function [i, u, F] = solveLinearisedWesterveltMultiLevelBoundaryExcitation(elements, omega, beta, gamma, x0, dx, nIterations, nHarmonics, linPointIsSolution)

n = size(elements.points,1);
N = nIterations;
h = zeros(n,1);
u = zeros(N,N,n);


F = zeros(N, n);
for i=1:N
    for j=0:min((i-1),nHarmonics-1)
        p_m_eta = zeros(1,n);
        p_m = zeros(1,n);
        % the first iteration has just the excitation on the right hand side
        % index 1 is the zero-th solution
        if i>1 && j>0
            % for the linearised operator the right hand side is a bit more
            % complex

            % ---------- First sum ----------
            % sum_{l=0}^j u_l * u_{j-l}
            for l = 0:j
                p_m_eta = p_m_eta + ...
                    squeeze(x0.u0(l+1,:)) .* ...
                    squeeze(u(i-1,(j-l)+1,:)).';
                 p_m = p_m + ...
                    squeeze(x0.u0(l+1,:)) .* ...
                    squeeze(x0.u0((j-l)+1,:));

            end

            % ---------- Second sum ----------
            % 2 * sum_{r=j}^{2*n - j} conj(u_r) * u_{r+j}
            for r = j:2:(2*(nHarmonics-1) - j)
                minusidx = (r-j)/2;
                plusidx = (r+j)/2;

                p_m_eta = p_m_eta + ...
                    conj(x0.u0(minusidx+1,:)) .* ...
                    squeeze(squeeze(u(i-1,plusidx+1,:)).') + ...
                    conj(squeeze(u(i-1,minusidx+1,:)).') .* ...
                    squeeze(x0.u0(plusidx+1,:));

                  p_m = p_m + 2 * ...
                    conj(squeeze(x0.u0(minusidx+1,:))) .* ...
                    squeeze(x0.u0(plusidx+1,:));
            end


        end
        cs = 1./(x0.s0.' + 1i.*j.*omega);
        cf = j.^2.*omega.^2.*cs;
        F(j+1,:) = -x0.eta0.'.*(p_m_eta).*cf; % et0 part
        F(j+1,:) = F(j+1, :) - dx.deta.'.*cf./2.*p_m; % deta part (we use the precomputed stuff from u0, this takes into account everything)
        F(j+1,:) = F(j+1, :) + dx.db.'.*cf.*x0.u0(j+1,:); %db part
        if (linPointIsSolution == false)
           F(j+1,:) = F(j+1, :) + dx.ds.'.*cs.*(x0.laplaceu0(j+1,:)); %ds part 
        else
           F(j+1,:) = F(j+1, :) + dx.ds.'.*cs.*(-x0.kappa0(:,j+1).'.*j^2.*x0.u0(j+1,:) - x0.F(j+1,:)); %ds part (here F(j+1,:) is correct)
        end
        F(j+1,elements.boundaryIdx) = 0;

        u(i,j+1,:) = solveHelmholtzCondensedC(elements, j*omega, gamma, j^2.*x0.kappa0(:,j+1), beta, F(j+1,:).', dx.excitation(:,j+1), n, elements.K, elements.rowK, elements.colK, elements.M_t, elements.tBM);

    end

end


end