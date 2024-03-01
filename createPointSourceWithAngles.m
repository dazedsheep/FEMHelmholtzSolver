function [source] = createPointSourceWithAngles(elements, sourceLocation, gamma, angles)
source = zeros(size(elements.points,1),1);

for j = 1:size(sourceLocation,2)
for i=1:size(elements.points,1)
    source(i,1) = source(i,1) + regularizedDirac(gamma,norm(elements.points(i,:).' - sourceLocation(:,j),2)).*angles(j);
end
end
end

