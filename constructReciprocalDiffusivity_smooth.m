function [b] = constructReciprocalDiffusivity(elements, diffusivityDomain, diffusivityPhantoms, centers, radii)
%% Squared speed of sound (space dependent)
% this does not yet have the refraction for phantoms included...

b = zeros(size(elements.points,1),1);

b(:) = 1./diffusivityDomain;

for j=1:size(centers,2)
    if(radii(j) == 0)
        % this is a point source
        % find nearest node to impose our point source
        [v,pcenterIdx] = min(sum((elements.points - centers(:,j)').^2,2)); 
        b(pcenterIdx) =  1./diffusivityPhantoms(j);
    else
        % this is a "disc" source
        if abs(diffusivityPhantoms(j)) > 0
            for i=1:size(elements.points,1)
                    b(i) = b(i) + (1./diffusivityPhantoms(j) -  1./diffusivityDomain).*regularizedDirac(radii(j), norm(elements.points(i,:) - centers(:,j)',2));
            end
        end
    end
end


end