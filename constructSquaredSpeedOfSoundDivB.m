function [s] = constructSquaredSpeedOfSoundDivB(elements, speed_of_sound, diffusivityDomain, diffusivityPhantoms, centers, radii)
%% Squared speed of sound (space dependent)
% this does not yet have the refraction for phantoms included...

s = zeros(size(elements.points,1),1);

s(:) = speed_of_sound^2./diffusivityDomain;

for j=1:size(centers,2)
    if(radii(j) == 0)
        % this is a point source
        % find nearest node to impose our point source
        [~,pcenterIdx] = min(sum((elements.points - centers(:,j)').^2,2)); 
        s(pcenterIdx) =  speed_of_sound.^2./diffusivityPhantoms(j);
    else
        % this is a "disc" source
        if abs(diffusivityPhantoms(j)) > 0
            for i=1:size(elements.points,1)
                if norm(elements.points(i,:) - centers(:,j)',2) < radii(j) 
                    s(i) = speed_of_sound.^2./diffusivityPhantoms(j);
                end
            end
        end
    end
end


end