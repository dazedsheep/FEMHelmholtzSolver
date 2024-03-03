function [source] = constructHalfPlanePointSource(elements, sourceLocation, gamma, directionOfPropagation)

source = createPointSource(elements, sourceLocation, gamma);

idx = abs(x) <= 2*gamma;

end
