function [elements] = initializeMultiLeveLSolver(meshSize, domain)
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

end