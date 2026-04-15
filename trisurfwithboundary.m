function trisurfwithboundary(elements, y)

figure, 
trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2), y, 'facecolor', 'interp'); 
shading interp;

hold on;

r = 0.2;
theta1 = pi/2;           % start angle (adjust as needed)
theta2 = 2*pi;       % end angle (adjust as needed)

theta = linspace(theta1, theta2, 200);

% Parametric equation of the arc
x_arc = r * cos(theta);
y_arc = r * sin(theta);

% Plot the boundary curve (slightly lifted in z for visibility)
z_arc = max(y(:)) * ones(size(x_arc));

plot3(x_arc, y_arc, z_arc, 'k', 'LineWidth', 2);

hold off;

end

