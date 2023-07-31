%% plain gradient descent for identifying sources and their values
m = 2;
sourcePoints = [];
sourceValues = [];

f = zeros(size(elements.points,1),1);
y = U(2,elements.bedges(:,1));
g = f;
h = f;
f(elements.bedges(:,1)) = -1i.*m.*kappa.*y;
% K*(-y) 
%te = solveHelmholtzVectorizedTmp(elements, m*w, m*kappa, -1/speed_of_sound, zeros(size(elements.points,1),1), f, g, size(elements.points,1));

%figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2), real(te), 'facecolor', 'interp'); shading interp;
%figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2), real(U(2,:)), 'facecolor', 'interp'); shading interp;


for j = 1:2
    
    te = solveHelmholtzVectorizedTmp(elements, m*w, m*kappa, -1/speed_of_sound, zeros(size(elements.points,1),1), f, g, size(elements.points,1));

    [v, pc] = max(abs(te));
    center = elements.points(pc,:)';
    sourcePoints = [sourcePoints center];
    % try to maximize the centers source value ||Kx - y||^2, where x are our
    % source values
    v = 1;
    F = constructF(elements, sourcePoints, zeros(size(sourcePoints,2),1), [sourceValues v]);
    v = [sourceValues v];
    fun = @(x) solveHelmholtzVectorizedTmp(elements, m*w, m*kappa, 1/speed_of_sound, constructF(elements, sourcePoints, zeros(size(sourcePoints,2),1), [x]).*m^2.*kappa^2.*coupling(m,:), g, g, size(elements.points,1));
    sigma = 10^-4;
    for k = 1:200
        lRate = 1000./min(abs(coupling(m,:)));
        % gradient
        h = 0.001;
        x = v;
        for i=1:size(x,2)
            e_i = zeros(1,size(x,2));
            e_i(i) = 1;
            y0 = fun(x-h*e_i);
            z = (y0(elements.bedges(:,1)).' - y);
            val0 = z*z';
    
    %         y0 = fun(x);
    %         z = (y0(elements.bedges(:,1)).' - y);
    %         val0 = z*z';
            y1 = fun(x+h*e_i);
            z = (y1(elements.bedges(:,1)).' - y);
            val1 = z*z';
            df(i,k) = 1/(2*h) *(val1 - val0);  
        end
        
        direction = -df(:,k);
        y_current = fun(x);
        loss = (y_current(elements.bedges(:,1)).' - y)*(y_current(elements.bedges(:,1)).' - y)';
        
        % line search / amijo conditions
        while 1
            y_new = fun(x + lRate*direction');
            loss_new = (y_new(elements.bedges(:,1)).' - y)*(y_new(elements.bedges(:,1)).' - y)';
            if (loss_new < loss + sigma*lRate*(direction'*direction))
                break;
            end
            lRate = lRate/10;
        end
    
        v = x - lRate*df(:,k)';
    end
    
    sourceValues = v;
    f(elements.bedges(:,1)) = 1i.*m.*kappa.*(y_new(elements.bedges(:,1)).' - y);


end
% uu = solveHelmholtzVectorizedTmp(elements, m*w, m*kappa, 1/speed_of_sound, F, g, g, size(elements.points,1));
% f(elements.bedges(:,1)) = -1i.*m.*kappa.*(uu(elements.bedges(:,1)).' - y);
% te2 = solveHelmholtzVectorizedTmp(elements, m*w, m*kappa, -1/speed_of_sound, zeros(size(elements.points,1),1), f, g, size(elements.points,1));
