function [elements] = createMesh(domain, H_max, H_min, H_edges)

model = createpde;
geom  = decsg(domain');
g = geometryFromEdges(model,geom);
%pdegplot(g,"VertexLabels","on","EdgeLabels","on")
generateMesh(model,'Hmax',H_max,'Hmin',H_min,'Hedge', {[1:g.NumEdges],H_edges}, 'GeometricOrder','Linear');
[p,e,t] = meshToPet(model.Mesh);
elements.tri = t(1 : 3, :)';
elements.points = p';
elements.edges = e([1 2 5],:)';
elements.T = triangulation(elements.tri(:,1 : 3),elements.points);
elements.otri = t;
elements.opoints = p;
%figure, pdegplot(model,'EdgeLabels','on');
%figure, pdeplot(model);
end