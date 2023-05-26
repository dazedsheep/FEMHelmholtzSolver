function [val] = h_func_ex2(x, y)

kappa = 20;
a = kappa.*cos(pi/2);
b = kappa.*sin(pi/2);
beta = 1;

bcenter = [1/2,1/2];
origu = @(x,y,a,b) (cos(a.*x+b.*y)+1i.*sin(a.*x + b.*y));
gradu = @(x, y, a, b) [-a.*sin(a.*x + b.*y) + 1i.*a.*cos(a.*x + b.*y), -b.*sin(a.*x + b.*y) + 1i.*b.*cos(a.*x + b.*y)];

nv = [x,y] - bcenter;
n_nv = nv./sqrt(sum(nv.^2,2));

val = diag(gradu(x,y,a,b)*n_nv') + 1i.*beta.*kappa.*origu(x,y,a,b);

end

