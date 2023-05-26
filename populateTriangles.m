function [triangles] = populateTriangles(elements)
    
for i=1:size(elements.tri,1)
    triangles{i}.triangleVertices = [elements.points(elements.tri(i,1),:); elements.points(elements.tri(i,2),:); elements.points(elements.tri(i,3),:)].';
end

end