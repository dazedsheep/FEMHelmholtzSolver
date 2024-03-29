function [i, boundaryIndices, elements, p, F, coupling] = solveForwardMultiLevelC(omega, gamma, beta, waveNumber, center, radii, sourceValues, domain, excitationPoints, excitationPointsSize, excitationPower, refractionIndex, speed_of_sound, diffusitivity, N, minHarmonics, threshold, massDensity, meshSize)

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
% not used
g = zeros(1,n);

% consider that 1/c^2 is already included
% select = zeros(n,1);
% select(200:1000) = 1;
% excitation = 10000*select; 

% 1-2*\eta*u > 0, otherwise the Westervelt equation degenerates
%excitation = 100*1/(2*max(sourceValues))*ones(n,1); %constant excitation?

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
                if sourceValues(j) > 0
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
        
        u(i,j,:) = solveHelmholtzCVectorizedKappaSampled(elements, (m+1)*omega, (m+1).*kappa, gamma,  beta, -F(j,:).', h, g, n);
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