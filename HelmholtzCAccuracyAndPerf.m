clear all
frequ = 1;
syms x y 
omega = 50;
kappa(x,y) = omega*((x+0.1)/(x+0.5));
origu(x,y) = (cos(kappa(x,y).*x+kappa(x,y).*y)+1i.*sin(kappa(x,y).*x + kappa(x,y).*y));

realf(x,y) = -laplacian(origu,[x,y]) - kappa(x,y)^2*origu(x,y);
gradu(x,y) = gradient(origu, [x,y]);

domainX = [0,1];
domainY = [0,1];

% --> Quadrangle format 1 x 10: [3 4 x_UpLeft x_UpR x_LoR x_LoLeft y_UpLeft y_UpR y_LoR y_LoLeft]
% --> Circle format 1 x 4: [1 x_c y_c radius]

% use the built-in meshing of MATLAB, we just have to cast the mesh into
% our structure then
H_max = 0.005;   
H_min = 0.005;
H_edges = 0.005;
domain = [3 4 0 1 1 0 1 1 0 0];
elements = createMesh(domain, H_max, H_min, H_edges);

figure, triplot(elements.T);

% use pdegplot to figure out the edge labels!!
%figure, pdegplot(model,'EdgeLabels','on');
elements.nr_edges = 1:4; 
elements.bedges = elements.edges(find(ismember(elements.edges(:,3),elements.nr_edges)),:);  % in our case these are all edges
elements.nodeIndex = elements.tri;

% populate triangles, this is still needed... (not nice)
elements.triangles = populateTriangles(elements);

beta = 1;
gamma = 0;
f = -realf;
innerpointsId = unique(elements.nodeIndex);
fI = double(f(elements.points(innerpointsId,1),elements.points(innerpointsId,2)));
kappaI = double(kappa(elements.points(innerpointsId,1),elements.points(innerpointsId,2)));
% CAREFUL!! - the h function has to be adapted according to the domain!
h = @(x,y) h_func_ex1c(x,y, beta, origu, gradu, omega, gamma);
for j=1:size(elements.bedges, 1)
     hI(elements.bedges(j,1)) = double(h(elements.points(elements.bedges(j,1), 1), elements.points(elements.bedges(j,1), 2)));
end
%
g = @(x,y) origu(x,y);
n = size(elements.points,1);
% tic
% U = solveHelmholtzC(elements, omega, kappa, gamma, beta, f, h, g, n);
% toc
tic
U = solveHelmholtzCVectorizedKappaSampled(elements, omega, kappaI, gamma, beta, fI, hI, g, n);
toc