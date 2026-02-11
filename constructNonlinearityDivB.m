function [f] = constructNonlinearityDivB(elements, massDensity, speed_of_sound, diffusivityDomain, diffusivityPhantoms, centers, radii, values, sourceValueDomain, convection)

if convection == true
    f  = (1+sourceValueDomain).*ones(size(elements.points,1),1)./(2.*massDensity.*speed_of_sound.*diffusivityDomain);
else
    f = zeros(size(elements.points,1),1);
end
% construct our f
for j=1:size(centers,2)
    if(radii(j) == 0)
        % this is a point source
        % find nearest node to impose our point source
        [v,pcenterIdx] = min(sum((elements.points - centers(:,j)').^2,2)); 
        f(pcenterIdx) =  (1+1./2.*values(j))./(2.*massDensity.*speed_of_sound.*diffusivityPhantoms(j));
    else
        % this is a "disc" source
        if abs(values(j)) > 0
            for i=1:size(elements.points,1)
                if norm(elements.points(i,:) - centers(:,j)',2) < radii(j) 
                    f(i) = (1+ 1./2.*values(j))./(2.*massDensity.*speed_of_sound.*diffusivityPhantoms(j));
                end
            end
        end
    end
end


end