function [i, u, F] = solveAdjointLinearisedWesterveltMultiLevelBoundaryExcitation(elements, omega, beta, gamma, x0, yobs, nIterations, nHarmonics)
% yobs has to be defined on the hole of \overline{\Omega}, just fill stuff
% with 0
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
u2_tt = zeros(N, n);
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
        % the boundary observation does not need to be scaled as we use the
        % normalised system!
        u(i,j+1,:) = solveHelmholtzCondensedC(elements, j*omega, gamma, conj(j^2.*x0.kappa0(:,j+1)), beta, F(j+1,:).',conj(cs).'.*yobs(j+1,:).', n, K, rowK, colK, M_t, tBM);
    end

end


end