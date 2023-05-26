function [elements] = domainTriangulation(domainX, domainY, n)
% This creates an uniform triangulation of a rectangular domain
% 
% domainX               ... array of lower and upper limit of X domain
% domainY               ... array of lower and upper limit of Y domain
% n                     ... number of nodes on axis

% raster on X and Y
rasterX = linspace(domainX(1), domainX(2),n);
rasterY = linspace(domainY(1), domainY(2),n);

% create a mesh
[X,Y] = meshgrid(rasterX,rasterY);

% gather the mesh points into a matrix
P = [X(:) Y(:)];

% define the dirichlet boundary node indices
boundaryNodes = [1:n,(n+1):n:(n*n-n),(2*n):n:(n*n-n),(n*n-n+1):1:(n*n)];

% triangles
DT = zeros(n*n,3);

ind = 1;

for i = 1:n-1
    for j = 1:n-1
        % two triangles meet at a node
        node         = (i-1)*n+j+1;              
        % triangle 1 (upper left)
        DT(ind,:)   = [node node-1 node+n];    
        % triangle 2 (lower right)
        DT(ind+1,:) = [node+n-1 node+n node-1];
        
        ind = ind+2;
    end
end
bnds = [boundaryNodes, 1];
edges = zeros(size(bnds,2)-1, 2);
for i=1:(size(bnds,2)-1)
    edges(i, 1) = bnds(i);
    edges(i, 2) = bnds(i+1);
end

% get the indices of the boundary nodes
% the triangle vertices are ordered counter clockwise in DT
k = 1;
for i=1:size(DT,1)

    elements.triangles{i}.triangleVertices = [P(DT(i,1),:); P(DT(i,2),:); P(DT(i,3),:)].';
    
    % are there nodes which form a boundary edge?
    isec = intersect(DT(i,:), boundaryNodes,'stable');
    % there have to be at least 2 nodes
    elements.triangles{i}.be = zeros(3,2);
    if (size(isec,2) > 1)
        
        elements.boundary{i}.bn  = isec;
        
        
        for j = 1:(size(isec,2)-1)
                
                tMI = inv(getTransformationMatrix(elements.triangles{i}.triangleVertices));
                if j==2
                    origTriVertices = tMI * ([P(isec(j+1),:)',P(isec(1),:)'] - elements.triangles{i}.triangleVertices(:,1));
                    x = P(isec(j+1),:);
                    y = P(isec(1),:);
                    idx = isec(j+1);
                    idy = isec(1);
                else
                    origTriVertices = tMI * ([P(isec(j),:)',P(isec(j+1),:)'] - elements.triangles{i}.triangleVertices(:,1));
                    x = P(isec(j),:);
                    y = P(isec(j+1),:);
                    idx = isec(j);
                    idy = isec(j+1);
                end
                v =  origTriVertices(:,2) - origTriVertices(:,1);
                
                if (abs(v(1)) > 0 && abs(v(2))>0)
                    continue
                end
                
                if (abs(v(1)) > 0)
                    elements.triangles{i}.be(2,:) = y - x;
                else
                    elements.triangles{i}.be(1,:) = y - x;
                end
                elements.edges(k,1) = idx;
                elements.edges(k,2) = idy;
                edges(k,1) = idx;
                edges(k,2) = idy;
                vec(k,:) = y - x;
                pnt(k,:) = x;
                elements.triangles{i}.bPoint(j,:) = x;
                k = k + 1;
                
        end
    end

end



elements.nodeIndex = DT;
elements.points = P;
elements.boundary = zeros(0,0);
elements.boundaryNeumann = boundaryNodes;
elements.tri = DT;
elements.bedges = edges;
% figure, quiver(pnt(:,1), pnt(:,2), vec(:,1), vec(:,2));
% 
% figure, plot(pnt);
% figure, triplot(DT, P(:,1), P(:,2));
% hold on
% plot(P(boundaryNodes(:), 1),P(boundaryNodes(:), 2),'*')
end

