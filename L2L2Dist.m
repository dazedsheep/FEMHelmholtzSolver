function [d] = L2L2Dist(elements, timeMesh, u, v)
diff = abs(u-v).^2;
[~, d] = integrate_fun_trimesh(elements.opoints, elements.otri, trapz(timeMesh.', diff, 1));
d = sqrt(d);
end

