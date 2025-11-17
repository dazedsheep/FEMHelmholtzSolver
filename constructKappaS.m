function  [kappa] = constructKappaS(elements, diffusivity, squared_speed_of_sound, omega, refractionIndex, sources, sourcesRadii, sourceValues, L, M, N)
%% Complex wave number
% with space dependent speed of sound and diffusvity

kappa = zeros(size(elements.points,1),L);
for j = 1:L
    [m,n] = linTo2D(j,M,N);
    kappa(:,j) = omega(m,n)./sqrt(squared_speed_of_sound + 1i*omega(m,n)*diffusivity(1));
end

for j=1:size(sources,2)
    if(sourcesRadii(j) == 0)
        % this is a point source
        % find nearest node to impose our point source
        [v,pcenterIdx] = min(sum((elements.points - sources(:,j)').^2,2)); 
        kappa(pcenterIdx,:) = omega./sqrt((squared_speed_of_sound(pcenterIdx)./refractionIndex(j)) + 1i.*omega.*diffusivity(j+1));
    else
        % this is a "disc" source
        for i=1:size(elements.points,1)
            if norm(elements.points(i,:) - sources(:,j)',2) < sourcesRadii(j) 
                if abs(sourceValues(j)) > 0
                    kappa(i,:) = omega(:)./sqrt((squared_speed_of_sound(i)./refractionIndex(j)) + 1i.*omega(:).*diffusivity(j+1));
                end
            end
        end
    end
end

end