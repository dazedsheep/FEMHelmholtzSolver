
% create the space dependent parameters
sourceValueDomain = 2; % B/A of domain

eta_values = [0.001]; % B/A of phantoms
eta_radii = [0.05];
eta_centers = [0.0;0.1];

s_values = [2.5]; % B/A of phantoms
s_radii = [0.05];
s_centers = [0.1;0.1];

b_values = [1.005]; % B/A of phantoms
b_radii = [0.05];
b_centers = [-0.1;-0.1];

eta = constructParameter(elements, eta_centers, eta_radii, eta_values,0);
s = constructParameter(elements, s_centers, s_radii, s_values,1);
b = constructParameter(elements, b_centers, b_radii, b_values,1);