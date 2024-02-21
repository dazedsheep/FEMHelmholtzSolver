function [source] = gaussianSource(elements, sourceLocation, sigma)
    source = zeros(size(elements.points,1),1);
    source(:,1) = 1/(sigma).*exp( -dot( (elements.points - repmat(sourceLocation.', size(elements.points,1),1)).', (elements.points - repmat(sourceLocation.', size(elements.points,1),1)).' )./(2.*sigma.^2));

end

