function [i, u, F] = solveLinearisedWesterveltMultiLevelBoundaryExcitation(elements, omega, beta, gamma, x0, dx, nIterations, nHarmonics, linPointIsSolution)

if (linPointIsSolution == false)
    error('not supported');
end

n = size(elements.points,1);
N = nIterations;
h = zeros(n,1);


u = zeros(N,N,n);

% extract nodes numbers of the 3 vertices of each triangle
n1x = elements.points(elements.tri(:,1),1).';
n1y = elements.points(elements.tri(:,1),2).';
n2x = elements.points(elements.tri(:,2),1).';
n2y = elements.points(elements.tri(:,2),2).';
n3x = elements.points(elements.tri(:,3),1).';
n3y = elements.points(elements.tri(:,3),2).';
m = size(elements.tri,1);
n1 = [n1x;n1y];
n2 = [n2x;n2y];
n3 = [n3x;n3y];

% compute the element wise transformation matrices
A = [n2(1,:) - n1(1,:), n3(1,:) - n1(1,:); n2(2,:) - n1(2,:), n3(2,:) - n1(2,:)];
B(1,1,:) = A(1,1:m);
B(1,2,:) = A(1,(m+1):2*m);
B(2,1,:) = A(2,1:m);
B(2,2,:) = A(2,(m+1):2*m);

% we need the inverse of each these matrices
C = pageinv(B);

% determinant of 3x3 matrix vectorized over all triangles
determinant = n1x .* n2y + n2x .* n3y  + n3x .* n1y - n1x .* n3y - n2x .* n1y - n3x .*n2y;
area = abs(determinant);

% calculate stiffness matrix
% the gradients of the P1 basis functions can be pre-computed
basisFuncGradients = [-1 -1; 1 0; 0 1];

% basisFuncGradients * inv(B(:,:,1)) * (basisFuncGradients * inv(B(:,:,1)))' * area(1)/2
K1 = pagemtimes(basisFuncGradients,C);

% finally, the stiffness matrix
K = pagemtimes(pagemtimes(K1,pagetranspose(K1)), reshape(area.*1/2, 1, 1, size(area,2)));

% calculate mass matrix
M_t = area/24 .* [2; 1; 1; 1; 2; 1; 1; 1; 2];

% local to global indices
rowK = elements.tri(:, [1 2 3 1 2 3 1 2 3]).';
colK = elements.tri(:, [1 1 1 2 2 2 3 3 3]).';

% compute length of boundary edges
e_Vec = elements.points(elements.bedges(:,1),:) - elements.points(elements.bedges(:,2),:);
e_len = sqrt(sum(e_Vec.^2,2));

% boundary element mass matrix
t_bM = e_len'/6 .* [2;1;1;2];

% local to global index for the boundary
brow = elements.bedges(:,[1 2 1 2]).';
bcol = elements.bedges(:,[1 1 2 2]).';

% sparse boundary mass matrix
tBM = sparse(brow, bcol, t_bM, size(elements.points,1),size(elements.points,1));

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
            % sum_{l=1}^j u_l * u_{j-l}
            for l = 1:j
                p_m_eta = p_m_eta + ...
                    squeeze(x0.u0(l+1,:)) .* ...
                    squeeze(u(i-1,(j-l)+1,:)).';
                 p_m = p_m + ...
                    squeeze(x0.u0(l+1,:)) .* ...
                    squeeze(x0.u0((j-l)+1,:));

            end

            % ---------- Second sum ----------
            % 2 * sum_{r=0}^{N-1-j} conj(u_r) * u_{r+j}
            for r = 0:((N-1)-j)
                p_m_eta = p_m_eta + ...
                    conj(x0.u0((r+j)+1,:)) .* ...
                    squeeze(squeeze(u(i-1,r+1,:)).') + ...
                    conj(squeeze(u(i-1,r+1,:)).') .* ...
                    squeeze(x0.u0((r+j)+1,:));
                  p_m = p_m + 2 * ...
                    conj(squeeze(x0.u0(r+1,:))) .* ...
                    squeeze(x0.u0((r+j)+1,:));
            end

        end

        
        F(j+1,:) = -x0.eta0.'.*(p_m_eta).*j^2.*x0.kappa0(:,j+1).'.*1./x0.b0.'; % et0 part
        F(j+1,:) = F(j+1, :) + dx.deta.'.*p_m; % deta part (we use the precomputed stuff from u0, this takes into account everything)
        F(j+1,:) = F(j+1, :) + dx.db.'.*j^2.*x0.kappa0(:,j+1).'.*1./x0.b0.'.*x0.u0(j+1,:); %db part
        F(j+1,:) = F(j+1, :) + dx.ds.'.*(1./(x0.s0.' + 1i.*j.*omega)).*(-x0.kappa0(:,j+1).'.*j^2.*x0.u0(j+1,:) + x0.F(j+1,:)); %ds part (here F(j+1,:) is correct)
        F(j+1,elements.boundaryIdx)  = 0;

        u(i,j+1,:) = solveHelmholtzCondensedC(elements, j*omega, gamma, j^2.*x0.kappa0(:,j+1), beta, F(j+1,:).', dx.excitation(:,j+1), n, K, rowK, colK, M_t, tBM);

    end

end


end