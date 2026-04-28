function trisurfwithboundary(elements, y, theta1, theta2, radius, color, linewidth, titleStr)

figure,

trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2), y, 'facecolor', 'interp'); 
shading interp;
title(titleStr);
hold on;
view(0,90)  
colorbar
xint = [0.14 0.22];
yint = [0.14 0.22];
annotation("textarrow",xint,yint,String="\Sigma");
set(gca,'fontname','Arial')  % Set it to times
set(gca, 'FontWeight', 'bold')

theta = linspace(theta1, theta2, 200);

% Parametric equation of the arc
x_arc = radius * cos(theta);
y_arc = radius * sin(theta);

% Plot the boundary curve (slightly lifted in z for visibility)
z_arc = max(y(:)) * ones(size(x_arc));

plot3(x_arc, y_arc, z_arc, color, 'LineWidth', linewidth);

hold off;

end

