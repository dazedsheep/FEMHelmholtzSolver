clear all
syms x y
domainX = [0,1];
domainY = [0,1];
omegas = 5:5:50;

edgelen = 0.005:0.005:0.1;

gridEdgeLen = 0.0005;
x = 0:gridEdgeLen:1;
y = x;

[X,Y] = meshgrid(x,y);

x1 = (0+gridEdgeLen):gridEdgeLen:(1-gridEdgeLen);
y1 = x1;
[X1,Y1] = meshgrid(x1,y1);

ctime = zeros(length(omegas), length(edgelen), 2);
cerror = zeros(length(omegas), length(edgelen), 3);


for o = 1:length(omegas)
    syms x y

    omega = omegas(o);
    kappa(x,y) = omega*((x+0.1)/(x+0.5));
    origu(x,y) = (cos(kappa(x,y).*x+kappa(x,y).*y)+1i.*sin(kappa(x,y).*x + kappa(x,y).*y));
    
    realf(x,y) = -laplacian(origu,[x,y]) - kappa(x,y)^2*origu(x,y);
    gradu(x,y) = gradient(origu, [x,y]);
    
    O = @(x,y) (cos(omega*((x+0.1)./(x+0.5)).*x+omega*((x+0.1)./(x+0.5)).*y)+1i.*sin(omega*((x+0.1)./(x+0.5)).*x + omega*((x+0.1)./(x+0.5)).*y));
    ru = O(X1,Y1);

% --> Quadrangle format 1 x 10: [3 4 x_UpLeft x_UpR x_LoR x_LoLeft y_UpLeft y_UpR y_LoR y_LoLeft]
% --> Circle format 1 x 4: [1 x_c y_c radius]
    for i= 1:length(edgelen)
        
        % use the built-in meshing of MATLAB, we just have to cast the mesh into
        % our structure then
        H_max = edgelen(i);   
        H_min = edgelen(i);
        H_edges = edgelen(i);
        domain = [3 4 0 1 1 0 1 1 0 0];
        elements = createMesh(domain, H_max, H_min, H_edges);
        
               
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
        tic
        U1 = solveHelmholtzCVectorized(elements, omega, kappa, gamma, beta, fI, hI, g, n);
        t1 = toc;
        
        tic
        U2 = solveHelmholtzCVectorizedKappaSampled(elements, omega, kappaI, gamma, beta, fI, hI, g, n);
        t2 = toc;
        
        U1xy = tri2grid(elements.opoints, elements.otri, U1, x1, y1 );
        U2xy = tri2grid(elements.opoints, elements.otri, U2, x1, y1 );

        ctime(o,i,1) = t1;
        ctime(o,i,2) = t2;
        cerror(o, i, 1) = sqrt(sum(sum( abs(( ru - U1xy )).^2.*gridEdgeLen^2)));
        cerror(o, i, 2) = sqrt(sum(sum( abs(( ru - U2xy )).^2.*gridEdgeLen^2)));
        cerror(o, i, 3) = sqrt(sum(sum( abs(( U1xy - U2xy )).^2.*gridEdgeLen^2)));

    end
end
%% plots
figure, plot(edgelen, squeeze(mean(ctime(:,:,1),1)), 'LineWidth', 1.5)
grid on
hold on
plot(edgelen, squeeze(mean(ctime(:,:,2),1)), 'LineWidth', 1.5, 'Color','red')
yscale log
xlabel("Mesh size (h)");
ylabel("Computation time [s]");
legend('Quadrature','Triangle constant \kappa');
set ( gca, 'xdir', 'reverse' )
%% kappa symbolic
figure, plot(edgelen, squeeze(cerror(:,:,1)), 'LineWidth', 1.5)
grid on
%hold on
%plot(edgelen, squeeze(mean(cerror(:,:,2),1)), 'LineWidth', 1.5, 'Color','red')
yscale log
xlabel("Mesh size (h)");
ylabel("L2 Error");
legendString = 'A = ' + string(omegas);
legend(legendString)
set ( gca, 'xdir', 'reverse' )

%%
figure, plot(edgelen, squeeze(cerror(:,:,2)), 'LineWidth', 1.5)
grid on
%hold on
%plot(edgelen, squeeze(mean(cerror(:,:,2),1)), 'LineWidth', 1.5, 'Color','red')
yscale log
xlabel("Mesh size (h)");
ylabel("L2 Error");
legendString = 'A = ' + string(omegas);
legend(legendString)
set ( gca, 'xdir', 'reverse' )

%% error between implementation variants
figure, plot(edgelen, squeeze(cerror(:,:,3)), 'LineWidth', 1.5)
grid on
%hold on
%plot(edgelen, squeeze(mean(cerror(:,:,2),1)), 'LineWidth', 1.5, 'Color','red')
yscale log
xlabel("Mesh size (h)");
ylabel("L2 Error");
legendString = 'A = ' + string(omegas);
legend(legendString)
set ( gca, 'xdir', 'reverse' )
