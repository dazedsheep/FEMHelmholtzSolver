function [i, boundaryIndices, elements, p, F, coupling] = solveForwardMultiLevelCVectorized(omega, gamma, beta, waveNumber, center, radii, sourceValues, domain, excitationPoints, excitationPointsSize, excitationPower, refractionIndex, speed_of_sound, diffusitivity, N, minHarmonics, threshold, massDensity, meshSize)

% domain triangulation
H_max = meshSize;   
H_min = meshSize;
H_edges = meshSize;
% our domain
bcenter = [domain(1),domain(2)];
brad = domain(3);

domain = [1, bcenter, brad];
elements = createMesh(domain, H_max, H_min, H_edges);

elements.nr_edges = 1:4; 
elements.bedges = elements.edges(find(ismember(elements.edges(:,3), elements.nr_edges)),:);  % in our case these are all edges
elements.nodeIndex = elements.tri;

% populate triangles, this is still needed... (not nice)
elements.triangles = populateTriangles(elements);

n = size(elements.points,1);

h = zeros(n,1);


excitation = getGridPoints(excitationPoints, excitationPointsSize, elements);
excitation = excitationPower * excitation;

f = zeros(n,1);
kappa = zeros(n,1);
kappa(:) = omega./sqrt(speed_of_sound.^2 + 1i*omega*diffusitivity);

% construct our f
for j=1:size(center,2)
    if(radii(j) == 0)
        % this is a point source
        % find nearest node to impose our point source
        [v,pcenterIdx] = min(sum((elements.points - center(:,j)').^2,2)); 
        f(pcenterIdx) = sourceValues(j)/(massDensity*speed_of_sound^2);
    else
        % this is a "disc" source
        for i=1:n
            if norm(elements.points(i,:) - center(:,j)',2) < radii(j) 
                f(i) = sourceValues(j)/(massDensity*speed_of_sound^2);
                if abs(sourceValues(j)) > 0
                    kappa(i) = omega/sqrt((speed_of_sound/refractionIndex(j))^2 + 1i*omega*diffusitivity);
                end
            end
        end
    end
end


waitbar_handle = waitbar(0,'Initializing waitbar...');

% now we iterate and apply the multilevel scheme that converges to the
% unique solution of the Westervelt equation

% we will solve N systems yielding N*(N-1)/2 helmholtz solutions but only
% the solution to the N-th systems is of interest
u = zeros(N,N,n);

% compute those elements apriori that do not change over the iterations

% initialize the needed data
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
t_bM = e_len'/6 .* [1;2;2;1]; 

% local to global index for the boundary
brow = elements.bedges(:,[1,2,1,2]);
bcol = elements.bedges(:,[1,1,2,2]);

% sparse boundary mass matrix
tBM = sparse(brow, bcol, t_bM, size(elements.points,1),size(elements.points,1));

waitbar(0, waitbar_handle, sprintf('%d of %d iterations done.', 0, N))
for i=1:N
   
    for j=1:i 
        m = j - 1;
        p_m = zeros(1,n);

        if i>1
            if j==2
                p_m = squeeze(2*sum(u(i-1,1:m,:).*u(i-1,(j-1):-1:(j-m),:),2)).';
            else
                p_m = squeeze(sum(u(i-1,1:m,:).*u(i-1,(j-1):-1:(j-m),:),2)).';
            end
        end

        if i > 1
            p_m = p_m + squeeze(sum(2*conj(u(i-1,(((j+2):2:(2*i-j))-j)/2,:)).*u(i-1,(((j+2):2:(2*i-j))+j)/2,:),2)).';
        end

        if j==1
            F(j,:) = excitation - 1/4.*f.*(m+1)^2.*kappa.^2.*p_m.';
        else
            F(j,:) = - 1/4.*f.*(m+1)^2.*kappa.^2.*p_m.';
        end

        coupling(j,:) = p_m;
        
        %u(i,j,:) = solveHelmholtzCVectorizedKappaSampled(elements, (m+1)*omega, (m+1).*kappa, gamma,  beta, -F(j,:).', h, g, n);
        u(i,j,:) = solveHelmholtzCondensedC(elements, (m+1)*omega, gamma, (m+1)*kappa, beta, -F(j,:).', h, n, K, rowK, colK, M_t, tBM);

    end
    waitbar(i/N, waitbar_handle, sprintf('%d of %d iterations done.', i, N))
    if i > minHarmonics
        if sum(abs(u(i-1,minHarmonics,:) - u(i,minHarmonics,:))) < threshold
            waitbar(i/N, waitbar_handle, sprintf('Threshold reached in the %d-th harmonic.', minHarmonics))
            break
        end
    end
end
p = u;
boundaryIndices = elements.bedges(:,1);

end