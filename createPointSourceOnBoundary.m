function [source] = createPointSourceOnBoundary(elements, sourceLocation, sourceSize)

source = zeros(size(elements.points,1),1);

% fetch the boundary points nearest to the source
% boundaryPointsIdx = zeros(size(sourceLocation,2),1);

for j = 1:size(sourceLocation,2)
    dists = repmat(sourceLocation(:,j)',size(elements.points(elements.bedges(:,1),:),1),1)  - elements.points(elements.bedges(:,1),:);

    myDistances = diag(dists*dists');
    % first find the nearest boundary point to the center of the source
    [~, idx] = min(myDistances);

    % distances on the boundary
    dists = repmat(elements.points(elements.bedges(idx,1),:), size(elements.points(elements.bedges(:,1),:),1),1) - elements.points(elements.bedges(:,1),:);
    myDistances = sqrt(diag(dists*dists.'));

    source(elements.bedges(:,1)) = source(elements.bedges(:,1)) + regularizedDirac(sourceSize, myDistances.').';

end

end

