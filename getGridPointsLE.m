function [gridPoints] = getGridPointsLE(elements, points, radii)
n = size(elements.points,1);
gridPoints = zeros(n,1);
% construct our f
for j=1:size(points,2)
    if(radii(j) == 0)
        % this is a point source
        % find nearest node to impose our point source
        [v,pcenterIdx] = min(sum((elements.points - points(:,j)').^2,2)); 
        gridPoints(pcenterIdx) = 1;
    else
        % this is a "disc" source
        for i=1:n
            if norm(elements.points(i,:) - points(:,j)',2) < radii(j) 
                gridPoints(i) = 1;
            end
        end
    end
end

end

