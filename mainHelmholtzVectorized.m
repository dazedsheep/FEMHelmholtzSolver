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
H_max = 0.005;   
H_min = 0.005;
H_edges = 0.005;
domain = [1, bcenter, brad];
elements = createMesh(domain, H_max, H_min, H_edges);

figure, triplot(elements.T);

% use pdegplot to figure out the edge labels!!
elements.nr_edges = 1:4; 
elements.bedges = elements.edges(find(ismember(elements.edges(:,3),elements.nr_edges)),:);  % in our case these are all edges
elements.nodeIndex = elements.tri;

% populate triangles, this is still needed... (not nice)
elements.triangles = populateTriangles(elements);

 kappa = 40;
    beta = 1;

f = @(x,y) realf(x,y,kappa*cos(pi/2),kappa*sin(pi/2), kappa);
h = @(x,y) h_func_ex2(x,y);
g = @(x,y) origu(x,y,kappa*cos(pi/2),kappa*sin(pi/2));

n = size(elements.points,1);
innerpointsId = unique(elements.nodeIndex);
fI = double(f(elements.points(innerpointsId,1),elements.points(innerpointsId,2)));
for i=1:size(elements.bedges, 1)
    hI(elements.bedges(i,1)) = h(elements.points(elements.bedges(i,1), 1), elements.points(elements.bedges(i,1), 2));
end

tic
U = solveHelmholtz(elements, kappa,0, beta, f, h, g, n);
t1 = toc;

disp(['FEM Solver for the Helmholtz equation took: ', num2str(t1),'s']);

tic
U = solveHelmholtzVectorized(elements, kappa, beta, fI, hI, g, n);
t2 = toc;

disp(['Vectorized FEM Solver for the Helmholtz equation took: ' , num2str(t2),'s']);
disp(['Speedup factor: ', num2str(t1/t2)]);

% tic
% U = solveHelmholtzVectorizedGPU(elements, kappa, 0, kappa, beta, fI, hI, n);
% t3 = toc;
% 
% disp(['GPU FEM Solver for the Helmholtz equation took: ', num2str(t3),'s']);

figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2),real(U), 'facecolor', 'interp'); shading interp;
title("Real part of the FEM solution using linear lagrange elements.")
xlabel('x');
ylabel('y');
figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2),imag(U), 'facecolor', 'interp'); shading interp;
title("Imag part of the FEM solution using linear lagrange elements.")
xlabel('x');
ylabel('y');


%% domain: rectangle in 2D - testi mesti
% this compares the non-vectorized FEM solver against a vectorized version
clear all

origu = @(x,y,a,b) (cos(a.*x+b.*y)+1i.*sin(a.*x + b.*y));
gradu = @(x, y, a, b, k) [-a*sin(a*x+b*y)+1i*a*cos(a*x + b*y), -b*sin(a*x + b*y) + 1i*b*cos(a*x + b*y)];
realf = @(x,y,a,b, kappa) -(a^2+b^2)*(cos(a*x + b*y) + 1i*sin(a*x+b*y)) + kappa^2*(cos(a*x+b*y)+1i*sin(a*x+b*y));

% our domain
domain = [3 4 0 1 1 0 1 1 0 0];
% use the built-in meshing of MATLAB, we just have to cast the mesh into
% our structure then
edgelen = 0.0005:0.0005:0.05;

gridEdgeLen = 0.0005;
x = 0:gridEdgeLen:1;
y = x;

[X,Y] = meshgrid(x,y);

x1 = (0+gridEdgeLen):gridEdgeLen:(1-gridEdgeLen);
y1 =x1;
[X1,Y1] = meshgrid(x1,y1);



kappas = 10:10:100;

ctime = zeros(length(kappas), length(edgelen), 2);
cerror = zeros(length(kappas), length(edgelen), 1);

for k = 1:length(kappas)
    kappa = kappas(k);
    beta = 1;

    ru = origu(X1,Y1,kappa*cos(pi/2),kappa*sin(pi/2));

    % figure, surface(x, y, real(Uxy)); shading flat;
    
    f = @(x,y) realf(x,y,kappa*cos(pi/2),kappa*sin(pi/2), kappa);
    h = @(x,y) h_func_ex1(x,y, kappa);
    g = @(x,y) origu(x,y,kappa*cos(pi/2),kappa*sin(pi/2));

    for i = 1:length(edgelen)
        
        H_max = edgelen(i);   
        H_min = edgelen(i);
        H_edges = edgelen(i);
        elements = createMesh(domain, H_max, H_min, H_edges);
        
        %figure, triplot(elements.T);
        
        % use pdegplot to figure out the edge labels!!
        elements.nr_edges = 1:4; 
        elements.bedges = elements.edges(find(ismember(elements.edges(:,3),elements.nr_edges)),:);  % in our case these are all edges
        elements.nodeIndex = elements.tri;
        
        % populate triangles, this is still needed... (not nice)
        elements.triangles = populateTriangles(elements);
        
        
        
        n = size(elements.points,1);
        
        innerpointsId = unique(elements.nodeIndex);
        fI = double(f(elements.points(innerpointsId,1),elements.points(innerpointsId,2)));
        for j=1:size(elements.bedges, 1)
            hI(elements.bedges(j,1)) = h(elements.points(elements.bedges(j,1), 1), elements.points(elements.bedges(j,1), 2));
        end
        
        tic
        U = solveHelmholtz(elements, kappa,0, beta, f, h, g, n);
        t1 = toc;
        
        %disp(['FEM Solver for the Helmholtz equation took: ', num2str(t1),'s']);
        
        tic
        U = solveHelmholtzVectorized(elements, kappa, beta, fI, hI, g, n);
        t2 = toc;
        
        ctime(k,i,1) = t1;
        ctime(k,i,2) = t2;
        
        % interpolate the result to the finer grid for evaluating the L2Error
        
        Uxy = tri2grid(elements.opoints, elements.otri, U, x1, y1 );
       
        cerror(k, i, 1) = sqrt(sum(sum( abs(( ru - Uxy )).^2.*gridEdgeLen^2)));
        
        %disp(['Vectorized FEM Solver for the Helmholtz equation took: ' , num2str(t2),'s']);
        %disp(['Speedup factor: ', num2str(t1/t2)]);
        
        
        
        % tic
        % U = solveHelmholtzVectorizedGPU(elements, kappa, 0, kappa, beta, fI, hI, n);
        % t3 = toc;
        % 
        % disp(['GPU FEM Solver for the Helmholtz equation took: ', num2str(t3),'s']);
        
        
        % figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2),real(U), 'facecolor', 'interp'); shading interp;
        % title("Real part of the FEM solution using linear lagrange elements.")
        % xlabel('x');
        % ylabel('y');
        % figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2),imag(U), 'facecolor', 'interp'); shading interp;
        % title("Imag part of the FEM solution using linear lagrange elements.")
        % xlabel('x');
        % ylabel('y');
    end
end
%% plots
figure, plot(edgelen, squeeze(mean(ctime,1)), 'LineWidth', 1.5)
grid on
yscale log
xlabel("Edge length (h)");
ylabel("Computation time [s]");
legend('Non-optimized solver','Optimized solver');
set ( gca, 'xdir', 'reverse' )

figure, plot(edgelen, cerror, 'LineWidth', 1.5);
xlabel("Edge length (h)");
ylabel("L2 Error");
yscale log
grid on
set ( gca, 'xdir', 'reverse' )
legendString = '\kappa = ' + string(kappas);
legend(legendString)
%% origu plots for different kappa
kappa = 50;
gridEdgeLen = 0.0005;
x = 0:gridEdgeLen:1;
y = x;

[X,Y] = meshgrid(x,y);
u = origu(X,Y,kappa,kappa);
figure, surf(X,Y, real(u)); shading flat; shading interp;
%title("Real part of exact sultion.")
xlabel('x');
ylabel('y');
view([-50,35]);
%az -50 EL 35
figure, surf(X,Y, imag(u)); shading flat; shading interp;
%title("Imag part of exact sultion.")
xlabel('x');
ylabel('y');
view([-50,35]);