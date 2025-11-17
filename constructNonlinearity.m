function [f] = constructNonlinearity(elements, massDensity, sources, sourcesRadii, sourceValues, sourceValueDomain, convection)

if convection == true
    f  = (1+sourceValueDomain).*ones(size(elements.points,1),1)./(massDensity);
else
    f = zeros(size(elements.points,1),1);
end
% construct our f
for j=1:size(sources,2)
    if(sourcesRadii(j) == 0)
        % this is a point source
        % find nearest node to impose our point source
        [v,pcenterIdx] = min(sum((elements.points - sources(:,j)').^2,2)); 
        f(pcenterIdx) =  (1+1./2.*sourceValues(j))./(massDensity);
    else
        % this is a "disc" source
        if abs(sourceValues(j)) > 0
            for i=1:size(elements.points,1)
                if norm(elements.points(i,:) - sources(:,j)',2) < sourcesRadii(j) 
                    f(i) = (1+ 1./2.*sourceValues(j))./(massDensity);
                end
            end
        end
    end
end


end