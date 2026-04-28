function [ut] = dt(u, timeMeshh)
% takes into account periodicity of u
ut = zeros(size(u));
ut = [ (u(end,:) - u(1,:))/timeMeshh; 
        (u(3:end,:) - u(1:end-2,:))/(2*timeMeshh); 
        (u(1,:) - u(end,:))/timeMeshh ];
end

