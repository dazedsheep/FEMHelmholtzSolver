function [s] = constructSquaredSpeedOfSound(elements, speed_of_sound)
%% Squared speed of sound (space dependent)
% this does not yet have the refraction for phantoms included...

s = zeros(size(elements.points,1),1);

s(:) = speed_of_sound^2;

end