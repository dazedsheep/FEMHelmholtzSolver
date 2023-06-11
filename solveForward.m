function [boundaryIndices, elements, p, F, coupling] = solveForward(c, omega, waveNumber, center, radii, sourceValues, domain, nHarmonics)
%% this function generates boundary data for given point sources in a ball 
H_max = 0.005;   
H_min = 0.005;

% our domain
bcenter = [domain(1),domain(2)];
brad = domain(3);

domain = [1, bcenter, brad];
elements = createMesh(domain, H_max, H_min);

elements.nr_edges = 1:4; 
elements.bedges = elements.edges(find(ismember(elements.edges(:,3),elements.nr_edges)),:);  % in our case these are all edges
elements.nodeIndex = elements.tri;

% populate triangles, this is still needed... (not nice)
elements.triangles = populateTriangles(elements);


% first we calculate p1 using a given excitation

kappa = waveNumber;
beta = 1/c;

n = size(elements.points,1);

h = zeros(n,1);

% not used
g = zeros(1,n);

% consider that 1/c^2 is already included
excitation = 2*ones(n,1); %constant excitation?

p1 = solveHelmholtzVectorizedTmp(elements, omega, kappa, beta, -excitation, h, g, n);

% figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2), real(p1), 'facecolor', 'interp'); shading interp;
% title("Real part of p_1(x).")
% xlabel('x');
% ylabel('y');

% figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2), imag(p1), 'facecolor', 'interp'); shading interp;
% title("Imag part of the FEM solution using linear lagrange elements.")
% xlabel('x');
% ylabel('y');


% now compute the solutions for the harmonic frequencies

p = zeros(nHarmonics+1, n);
p(1,:) = p1;
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

figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2), f, 'facecolor', 'interp'); shading interp;
% title("Scatterers.")
F(1,:) = p1;
for j = 1:(nHarmonics)
    m = j + 1;
    % not nice
    p_i = zeros(1,n);
    % TODO: vectorize
     for i = 1:j
         p_i = p_i + p(i,:).*p(m-i,:);
     end
    F(m,:) = -f.*m^2.*kappa^2.*p_i';
    coupling(j,:) = p_i;
    p(m,:) = solveHelmholtzVectorizedTmp(elements, m*omega, m*kappa, beta, -f.*m^2.*kappa^2.*p_i', h, g, n);
end

%%
boundaryIndices = elements.bedges(:,1);
end