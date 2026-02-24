function [val] = innerProd(elements, timeMeshh, a, b)

[vals, ~] = integrate_fun_trimesh(elements.opoints, elements.otri, a.s .* b.s);
[valb, ~] = integrate_fun_trimesh(elements.opoints, elements.otri, a.b .* b.b);
[valeta, ~] = integrate_fun_trimesh(elements.opoints, elements.otri, a.eta .* b.eta);

% TODO inner product on X_u_j (H^1(0,T_j,H^k(\Omega)))
valu = 0;
for k = 1:size(a.u,1)
    [v, ~] = integrate_fun_trimesh(elements.opoints, elements.otri, a.u(k,:) .* b.u(k,:));
    valu = valu + timeMeshh*v; % weighted sum over time
end

val = vals + valb  + valeta + valu;

end

