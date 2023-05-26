function [elements] = createMesh(domain, H_max, H_min)

model = createpde;
geom  = decsg(domain');
geometryFromEdges(model,geom);
generateMesh(model,'Hmax',H_max,'Hmin',H_min,'GeometricOrder','Linear');
[p,e,t] = meshToPet(model.Mesh);
elements.tri = t(1 : 3, :)';
elements.points = p';
elements.edges = e([1 2 5],:)';
elements.T = triangulation(elements.tri(:,1 : 3),elements.points);
%figure, pdegplot(model,'EdgeLabels','on');
end

