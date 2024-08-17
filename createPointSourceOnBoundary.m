function [source] = createPointSourceOnBoundary(elements, sourceLocation, sourceSize, gamma)

source = zeros(size(elements.points,1),1);

% fetch the boundary points nearest to the source
% boundaryPointsIdx = zeros(size(sourceLocation,2),1);

for j = 1:size(sourceLocation,2)
    dists = repmat(sourceLocation(:,j)',size(elements.points(elements.bedges(:,1),:),1),1)  - elements.points(elements.bedges(:,1),:);
    
    myDistances = diag(dists*dists');
    % first find the nearest boundary point to the center of the source
    [dist, idx] = min(myDistances);
    
    % distances on the boundary
    dists = repmat(elements.points(elements.bedges(idx,1),:), size(elements.points(elements.bedges(:,1),:),1),1) - elements.points(elements.bedges(:,1),:);
    myDistances = sqrt(diag(dists*dists.'));
    
    source(elements.bedges(:,1)) = source(elements.bedges(:,1)) + regularizedDirac(gamma, myDistances.').';


end

%source(boundaryPointsIdx,1) = ones(size(boundaryPointsIdx,1),1);

% for j = 1:size(sourceLocation,2)
%     for i=1:size(elements.points,1)
%         %source(i,1) = source(i,1) + regularizedDirac(gamma,norm(elements.points(i,:).' - sourceLocation(:,j),2));
% 
%     end
% end

end

