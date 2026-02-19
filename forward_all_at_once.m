function F = forward_all_at_once(u, s, b, eta, M, K, Mb, Bx, By, Nx, Ny, omega, gamma, Sigma, boundaryNodes)
%
%
% Output structure:
%   F{j}.int   : interior residual (volume)
%   F{j}.bdry  : boundary residual (Robin)
%   F{j}.obs   : observation residual
%


nExp = 1;

F = cell(nExp,1);

for j = 1:nExp
    
    
    [nNodes, nHarm] = size(u);
    
    F_int  = zeros(nNodes,nHarm);
    F_bdry = zeros(length(boundaryNodes),nHarm);
    F_obs  = zeros(length(Sigma),nHarm);
    
    % harmonic convolution for u^2
    u_sq = zeros(nNodes,nHarm);
    for k = 0:nHarm-1
        for l = 0:k
            u_sq(:,k+1) = u_sq(:,k+1) ...
                + u(:,l+1).*u(:,k-l+1);
        end
    end
    
    for k = 0:nHarm-1
        
        uk = u(:,k+1);
        
        % -------- INTERIOR --------
        
        nonlinear = -(k*omega)^2 * ...
            ( M * ( b.*uk ...
            - eta.*u_sq(:,k+1) ) );
        
        diffusion = s .* (K*uk);
        
        damping   = 1i*k*omega * (K*uk);
        
        F_int(:,k+1) = nonlinear ...
                     + diffusion ...
                     + damping;
        
        
        % -------- BOUNDARY (Robin) --------
        
        robin_mass = gamma * (Mb * uk);
        
        ukx = Bx * uk;
        uky = By * uk;
     
        normal_flux = (ukx(boundaryNodes) .* Nx + uky(boundaryNodes) .* Ny); 
        F_bdry(:,k+1) = robin_mass(boundaryNodes) + normal_flux;
        
        
        % -------- OBSERVATION --------
        
        F_obs(:,k+1) = uk(Sigma);
        
    end
    
    F{j}.int  = F_int;
    F{j}.bdry = F_bdry;
    F{j}.obs  = F_obs;
end

end
