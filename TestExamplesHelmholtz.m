%%

clear all
origu = @(x,y,a,b) (cos(a.*x+b.*y)+1i.*sin(a.*x + b.*y));
gradu = @(x, y, a, b, k) [-a*sin(a*x+b*y)+1i*a*cos(a*x + b*y), -b*sin(a*x + b*y) + 1i*b*cos(a*x + b*y)];
realf = @(x,y,a,b, kappa) -(a.^2+b.^2).*(cos(a.*x + b.*y) + 1i.*sin(a.*x+b.*y)) + kappa.^2.*(cos(a.*x+b.*y)+1i.*sin(a.*x+b.*y));

kappa= 20;
beta = 1;
domainX = [0,1];
domainY = [0,1];
n = 128;
elements = domainTriangulation(domainX, domainY, n);
n = size(elements.points,1);

f = @(x,y) realf(x,y,kappa*cos(pi/2),kappa*sin(pi/2), kappa);
h = @(x,y) h_func_ex1(x,y);
g = @(x,y) origu(x,y,kappa*cos(pi/2),kappa*sin(pi/2));


% plot the numerical solution
U = solveHelmholtz(elements, kappa, beta, f, h, g, n);
figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2),real(U), 'facecolor', 'interp'); shading interp;
title("Real part of the FEM solution using linear lagrange elements.")
xlabel('x');
ylabel('y');

figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2),imag(U), 'facecolor', 'interp'); shading interp;
title("Imag part of the FEM solution using linear lagrange elements.")
xlabel('x');
ylabel('y');


[X,Y] = meshgrid(linspace(domainX(1),domainX(2)), linspace(domainY(1),domainY(2)));
u = origu(X,Y,kappa*cos(pi/2),kappa*sin(pi/2));
figure, surf(X,Y, real(u)); shading flat; shading interp;
title("Real part of exact sultion.")
xlabel('x');
ylabel('y');

figure, surf(X,Y, imag(u)); shading flat; shading interp;
title("Imag part of exact sultion.")
xlabel('x');
ylabel('y');

%% domain: rectangle [0,1]x[0,1] using built-in mesh generation
clear all

origu = @(x,y,a,b) (cos(a.*x+b.*y)+1i.*sin(a.*x + b.*y));
gradu = @(x, y, a, b, k) [-a*sin(a*x+b*y)+1i*a*cos(a*x + b*y), -b*sin(a*x + b*y) + 1i*b*cos(a*x + b*y)];
realf = @(x,y,a,b, kappa) -(a^2+b^2)*(cos(a*x + b*y) + 1i*sin(a*x+b*y)) + kappa^2*(cos(a*x+b*y)+1i*sin(a*x+b*y));

% --> Quadrangle format 1 x 10: [3 4 x_UpLeft x_UpR x_LoR x_LoLeft y_UpLeft y_UpR y_LoR y_LoLeft]
% --> Circle format 1 x 4: [1 x_c y_c radius]

% use the built-in meshing of MATLAB, we just have to cast the mesh into
% our structure then
H_max = 0.015;   
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

kappa= 20;
beta = 1;

f = @(x,y) realf(x,y,kappa*cos(pi/2),kappa*sin(pi/2), kappa);

% CAREFUL!! - the h function has to be adapted according to the domain!
h = @(x,y) h_func_ex1(x,y);
%
g = @(x,y) origu(x,y,kappa*cos(pi/2),kappa*sin(pi/2));
n = size(elements.points,1);
U = solveHelmholtz(elements, kappa, beta, f, h, g, n);

figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2),real(U), 'facecolor', 'interp'); shading interp;
title("Real part of the FEM solution using linear lagrange elements.")
xlabel('x');
ylabel('y');

