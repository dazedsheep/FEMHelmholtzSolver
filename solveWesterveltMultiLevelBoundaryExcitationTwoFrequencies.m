function [i, u, F] = solveWesterveltMultiLevelBoundaryExcitationTwoFrequencies(elements, squared_speed_of_sound, omega, beta, gamma, kappa, excitation, f, M, N, threshold)
%% Input formats
% elements: finite element structure
% omega: 2*pi*f (an array with all the frequencies including harmonics)
% gamma, beta: impedance boundary parameters (space dependent, sampled at
% the boundary points in "elements")
% kappa: space dependent wavenumbers for all considered combinations of fundamental frequencies (sampled on points in
% elements), format is L = (M+1)*(N+1) - 1 (harmonics) x elements
% excitation: boundary excitation(s) of the format elements x
% frequencies (indexed one dimensionally)
% f: space dependent nonlineartity parameter 
% nHarmonics: how many harmonics are we going to compute (total number) =
% M*N

n = size(elements.points,1);
L = M*N; % number of solution we have to compute in total

waitbar_handle = waitbar(0,'Initializing waitbar...');

% pre-alloc the matrix for the Helmholtz solutions; for now we live with
% the memeory it eats because we want to inspect the solutions in every
% iteration
u = zeros(L,L,n);

% prepare as much as possible for the conforming FEM Helmholtz solver

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

handleWaitbar = waitbar(0, waitbar_handle, sprintf('%d of %d iterations done.', 0, L));
F = zeros(L, n);
elapsedTime = -1;

% skip the first solution as it is anyway zero!
for i=2:L 
    % tell the user what harmonic we are computing
    
    sprintf('Computing (%d,%d) harmonics...\n', linTo2D(i,M,N));
    
    for j=1:i
        if i==2 && j==2
            tic
        end

        m = j - 1;
        p_m = zeros(1,n);

        % the first iteration has just the excitation on the right hand side
        if i>1
            p_m = squeeze(sum(u(i-1,1:m,:).*u(i-1,(j-1):-1:(j-m),:),2)).';
            p_m = p_m + 2.*squeeze(sum(conj(u(i-1,(((j+2):2:(2*i-j))-j)/2,:)).*u(i-1,(((j+2):2:(2*i-j))+j)/2,:),2)).';
        end

        F(j,:) = kappa(:,j).^2.*(-1/2 .* f./squared_speed_of_sound .* p_m.');

        u(i,j,:) = solveHelmholtzCondensedC(elements, omega(j), gamma, kappa(:,j), beta, F(j,:).', -excitation(:,j), n, K, rowK, colK, M_t, tBM);
     
        if i==2 && j==2
            elapsedTime = toc;
        end
    end
    
    waitbar(i/L, waitbar_handle, sprintf('%d of %d iterations done (est. time left %f s).', i, L, elapsedTime*(L*L - i*j)))
    %% TODO: fix threshold
    % if i > minHarmonics
    %     if sum(abs(u(i-1,minHarmonics,:) - u(i,minHarmonics,:))) < threshold
    %         waitbar(i/L, waitbar_handle, sprintf('Threshold reached in the %d-th harmonic.', minHarmonics))
    %         break
    %     end
    % end
end

close(handleWaitbar);

end