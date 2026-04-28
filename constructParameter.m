function  [param] = constructParameter(elements, centers, radii, values, minvalue)

param = ones(size(elements.points,1),1)*minvalue;

for j=1:size(centers,2)
    if(radii(j) == 0)
        % this is a point source
        % find nearest node to impose our point source
        [~,pcenterIdx] = min(sum((elements.points - centers(:,j)').^2,2)); 
        param(pcenterIdx) =  values(j);
    else
        % this is a "disc" source
        if abs(values(j)) > 0
            for i=1:size(elements.points,1)
                if norm(elements.points(i,:) - centers(:,j)',2) < radii(j) 
                    param(i) = values(j);
                end
            end
        end
    end
end

end