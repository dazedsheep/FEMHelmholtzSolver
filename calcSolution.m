function [u] = calcSolution(elements, timeMesh, U, omega)

% U includes all harmonics
time_func = exp(1i .* (1:size(U,1)).'.* omega.* timeMesh);

for i=1:size(U,1)
    uHarmonics(i, :, :) = time_func(i,:).'.* U(i,:);
end

u = real(squeeze(sum(uHarmonics,1)));
end