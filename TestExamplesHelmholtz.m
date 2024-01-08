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

domainX = [0,1];
domainY = [0,1];

% --> Quadrangle format 1 x 10: [3 4 x_UpLeft x_UpR x_LoR x_LoLeft y_UpLeft y_UpR y_LoR y_LoLeft]
% --> Circle format 1 x 4: [1 x_c y_c radius]

% use the built-in meshing of MATLAB, we just have to cast the mesh into
% our structure then
H_max = 0.1;   
H_min = 0.1;
H_edges = 0.1;
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

kappa = 2;
beta = 1;

f = @(x,y) realf(x,y,kappa,kappa, kappa);

% CAREFUL!! - the h function has to be adapted according to the domain!
h = @(x,y) h_func_ex1(x,y, kappa);
%
g = @(x,y) origu(x,y,kappa,kappa);
n = size(elements.points,1);
U = solveHelmholtz(elements, kappa, 0, beta, f, h, g, n);
% fI = f(elements.points(elements.nodeIndex(:,1),1),elements.points(elements.nodeIndex(:,1),2));
% for i=1:size(elements.bedges, 1)
%     hI(elements.bedges(i,1)) = h(elements.points(elements.bedges(i,1), 1), elements.points(elements.bedges(i,1), 2));
% end
% 
% U = solveHelmholtzVectorizedGPU(elements, kappa, 0, kappa, beta, fI, hI, n);
% 

figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2),real(U), 'facecolor', 'interp'); shading interp;
title("Real part of the FEM solution using linear lagrange elements.")
xlabel('x');
ylabel('y');

[X,Y] = meshgrid(linspace(domainX(1),domainX(2)), linspace(domainY(1),domainY(2)));
u = origu(X,Y,kappa,kappa);
figure, surf(X,Y, real(u)); shading flat; shading interp;
title("Real part of exact sultion.")
xlabel('x');
ylabel('y');

figure, surf(X,Y, imag(u)); shading flat; shading interp;
title("Imag part of exact sultion.")
xlabel('x');
ylabel('y');

%% Example with space dependent kappa>0 a.e. (or speed of propagation)
clear all
frequ = 1;
syms x y 
omega = 4;
kappa(x,y) = 4*((x)/(x+1));
origu(x,y) = (cos(kappa(x,y).*x+kappa(x,y).*y)+1i.*sin(kappa(x,y).*x + kappa(x,y).*y));

realf(x,y) = -laplacian(origu,[x,y]) - kappa(x,y)^2*origu(x,y);
gradu(x,y) = gradient(origu, [x,y]);

domainX = [0,1];
domainY = [0,1];

% --> Quadrangle format 1 x 10: [3 4 x_UpLeft x_UpR x_LoR x_LoLeft y_UpLeft y_UpR y_LoR y_LoLeft]
% --> Circle format 1 x 4: [1 x_c y_c radius]

% use the built-in meshing of MATLAB, we just have to cast the mesh into
% our structure then
H_max = 0.05;   
H_min = 0.05;
H_edges = 0.05;
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
U = solveHelmholtzCVectorized(elements, omega, kappa, gamma, beta, fI, hI, g, n);
toc
figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2),real(U), 'facecolor', 'interp'); shading interp;
title("Real part of the FEM solution using linear lagrange elements.")
xlabel('x');
ylabel('y');

[X,Y] = meshgrid(linspace(domainX(1),domainX(2)), linspace(domainY(1),domainY(2)));
%O = @(x,y) (cos(kappa(x,y).*x+kappa(x,y).*y)+1i.*sin(kappa(x,y).*x + kappa(x,y).*y));
O = @(x,y) (cos(4.*((x)./(x+1)).*x+4.*((x)./(x+1)).*y)+1i.*sin(4.*((x)./(x+1)).*x + 4.*((x)./(x+1)).*y));
u = O(X,Y);
figure, surf(X,Y, real(u)); shading flat; shading interp;
title("Real part of exact sultion.")
xlabel('x');
ylabel('y');

figure, surf(X,Y, imag(u)); shading flat; shading interp;
title("Imag part of exact sultion.")
xlabel('x');
ylabel('y');

%% Example with space dependent kappa>0 a.e. (or speed of propagation) - with sampled kappa for comparison
clear all
frequ = 1;
syms x y 
omega = 25;
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
% tic
% U = solveHelmholtzCVectorized(elements, omega, kappa, gamma, beta, fI, hI, g, n);
% toc
tic
U = solveHelmholtzCVectorizedKappaSampled(elements, omega, kappaI, gamma, beta, fI, hI, g, n);
toc
figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2),real(U), 'facecolor', 'interp'); shading interp;
title("Real part of the FEM solution using linear lagrange elements.")
xlabel('x');
ylabel('y');

[X,Y] = meshgrid(linspace(domainX(1),domainX(2)), linspace(domainY(1),domainY(2)));
%O = @(x,y) (cos(kappa(x,y).*x+kappa(x,y).*y)+1i.*sin(kappa(x,y).*x + kappa(x,y).*y));
O = @(x,y) (cos(omega*((x+0.1)./(x+0.5)).*x+omega*((x+0.1)./(x+0.5)).*y)+1i.*sin(omega*((x+0.1)./(x+0.5)).*x + omega*((x+0.1)./(x+0.5)).*y));
u = O(X,Y);
figure, surf(X,Y, real(u)); shading flat; shading interp;
title("Real part of exact sultion.")
xlabel('x');
ylabel('y');

figure, surf(X,Y, imag(u)); shading flat; shading interp;
title("Imag part of exact sultion.")
xlabel('x');
ylabel('y');

%%
omega = 50;

domainX = [0,1];
domainY = [0,1];
[X,Y] = meshgrid(linspace(domainX(1),domainX(2)), linspace(domainY(1),domainY(2)));

O = @(x,y) (cos(omega*((x+0.1)./(x+0.5)).*x+omega*((x+0.1)./(x+0.5)).*y)+1i.*sin(omega*((x+0.1)./(x+0.5)).*x + omega*((x+0.1)./(x+0.5)).*y));
u = O(X,Y);
figure, surf(X,Y, real(u)); shading flat; shading interp;
%title("Real part of exact sultion.")
view([-50,35]);
xlabel('x');
ylabel('y');
figure, surf(X,Y, imag(u)); shading flat; shading interp;
%title("Imag part of exact sultion.")
view([-50,35]);
xlabel('x');
ylabel('y');