%% domain: ball in 2D
% this compares the non-vectorized FEM solver against a vectorized version
clear all

origu = @(x,y,a,b) (cos(a.*x+b.*y)+1i.*sin(a.*x + b.*y));
gradu = @(x, y, a, b, k) [-a*sin(a*x+b*y)+1i*a*cos(a*x + b*y), -b*sin(a*x + b*y) + 1i*b*cos(a*x + b*y)];
realf = @(x,y,a,b, kappa) -(a^2+b^2)*(cos(a*x + b*y) + 1i*sin(a*x+b*y)) + kappa^2*(cos(a*x+b*y)+1i*sin(a*x+b*y));

% our domain
bcenter = [1/2,1/2];
brad = 1/2;

% use the built-in meshing of MATLAB, we just have to cast the mesh into
% our structure then
H_max = 0.015;   
H_min = 0.002;
domain = [1, bcenter, brad];
elements = createMesh(domain, H_max, H_min);

figure, triplot(elements.T);

% use pdegplot to figure out the edge labels!!
elements.nr_edges = 1:4; 
elements.bedges = elements.edges(find(ismember(elements.edges(:,3),elements.nr_edges)),:);  % in our case these are all edges
elements.nodeIndex = elements.tri;

% populate triangles, this is still needed... (not nice)
elements.triangles = populateTriangles(elements);

kappa = 20;
beta = 1;

f = @(x,y) realf(x,y,kappa*cos(pi/2),kappa*sin(pi/2), kappa);
h = @(x,y) h_func_ex2(x,y);
g = @(x,y) origu(x,y,kappa*cos(pi/2),kappa*sin(pi/2));

n = size(elements.points,1);

tic
U = solveHelmholtz(elements, kappa, beta, f, h, g, n);
t1 = toc;

disp(['FEM Solver for the Helmholtz equation took: ', num2str(t1),'s']);

tic
U = solveHelmholtzVectorized(elements, kappa, beta, f, h, g, n);
t2 = toc;

disp(['Vectorized FEM Solver for the Helmholtz equation took: ' , num2str(t2),'s']);
disp(['Speedup factor: ', num2str(t1/t2)]);

figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2),real(U), 'facecolor', 'interp'); shading interp;
title("Real part of the FEM solution using linear lagrange elements.")
xlabel('x');
ylabel('y');
figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2),imag(U), 'facecolor', 'interp'); shading interp;
title("Imag part of the FEM solution using linear lagrange elements.")
xlabel('x');
ylabel('y');