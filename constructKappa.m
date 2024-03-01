function  [kappa] = constructKappa(elements, diffusivity, speed_of_sound, omega, refractionIndex, sources, sourcesRadii, sourceValues, nHarmonics)

kappa = zeros(size(elements.points,1),nHarmonics);
for j = 1:nHarmonics
    kappa(:,j) = omega./sqrt(speed_of_sound.^2 + 1i*j*omega*diffusivity);
end
for j=1:size(sources,2)
    if(sourcesRadii(j) == 0)
        % this is a point source
        % find nearest node to impose our point source
        [v,pcenterIdx] = min(sum((elements.points - sources(:,j)').^2,2)); 
        kappa(pcenterIdx,:) = omega./sqrt((speed_of_sound./refractionIndex(j))^2+ 1i.*(1:nHarmonics).*omega*diffusivity);
    else
        % this is a "disc" source
        for i=1:size(elements.points,1)
            if norm(elements.points(i,:) - sources(:,j)',2) < sourcesRadii(j) 
                if abs(sourceValues(j)) > 0
                    kappa(i,:) = omega./sqrt((speed_of_sound/refractionIndex(j))^2 + 1i.*(1:nHarmonics).*omega*diffusivity);
                end
            end
        end
    end
end

end