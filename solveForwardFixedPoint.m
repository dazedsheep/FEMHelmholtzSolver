function [boundaryIndices, elements, p, F, coupling] = solveForwardFixedPoint(c, omega, gamma, waveNumber, center, radii, sourceValues, domain, N)


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
excitation = 1/max(sourceValues).*ones(n,1); %constant excitation?

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

% we will do a fixed point iteration up to the second harmonic
%
u = zeros(N,2,n);
coupling = zeros(N,2,n);

% i ... iteration
for i = 2:N
    m = 1;
    p_1 = squeeze(conj(u(i-1,1,:)).*u(i-1,2,:)).';
    f1 = excitation - 1/2.*f.*m^2.*kappa^2.*p_1.';
    coupling(i,m,:) = p_1;
    % u_1 in this iteration
    u(i,1,:) = solveHelmholtzVectorizedTmp(elements, m*omega, gamma, m*kappa, beta, -f1, h, g, n);
    
    
    m = 2;
    p_2 = squeeze(u(i,1,:).*u(i,1,:)).';
    f2 = -1/4.*f.*m^2.*kappa^2.*p_2.';
    coupling(i,m,:) = p_2;
    % update u2
    u(i,2,:) = solveHelmholtzVectorizedTmp(elements, m*omega, gamma, m*kappa, beta, -f2, h, g, n);
end

F = f1;
p = u;
boundaryIndices = elements.bedges(:,1);



end

