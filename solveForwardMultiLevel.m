function [boundaryIndices, elements, p, F, coupling] = solveForwardMultiLevel(c, omega, gamma, waveNumber, center, radii, sourceValues, domain, excitationPoints, excitationPointsSize, excitationPower, N)

% domain triangulation
H_max = 0.005;   
H_min = 0.005;
H_edges = 0.005;
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

kappa = waveNumber;
beta = 1/c;
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
% construct our f
for j=1:size(center,2)
    if(radii(j) == 0)
        % this is a point source
        % find nearest node to impose our point source
        [v,pcenterIdx] = min(sum((elements.points - center(:,j)').^2,2)); 
        f(pcenterIdx) = sourceValues(j);
    else
        % this is a "disc" source
        for i=1:n
            if norm(elements.points(i,:) - center(:,j)',2) < radii(j) 
                f(i) = sourceValues(j);
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
                for s=1:m
                    p_m = p_m + 2*squeeze(u(i-1,s,:).*u(i-1,j-s,:)).';
                end
            else
                for s=1:m
                    p_m = p_m + squeeze(u(i-1,s,:).*u(i-1,j-s,:)).';
                end
            end
        end

        if i > 1
            for s = (j+2):2:(2*i-j)
                p_m = p_m + squeeze(2*conj(u(i-1,(s-j)/2,:)).*u(i-1,(s+j)/2,:)).';
            end
        end
        if j==1
            F(j,:) = excitation - 1/4.*f.*(m+1)^2.*kappa^2.*p_m.';
        else
            F(j,:) = - 1/4.*f.*(m+1)^2.*kappa^2.*p_m.';
        end
        coupling(j,:) = p_m;
        u(i,j,:) = solveHelmholtzVectorizedTmp(elements, (m+1)*omega, gamma, (m+1)*kappa, beta, -F(j,:), h, g, n);

    end
    waitbar(i/N, waitbar_handle, sprintf('%d of %d iterations done.', i, N))
end
p = u;
boundaryIndices = elements.bedges(:,1);

end