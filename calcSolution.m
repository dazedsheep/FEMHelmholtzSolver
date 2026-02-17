function [u] = calcSolution(elements, timeMesh, U, omega)
    u = real( exp(1i*(0:size(U,1)-1).' * omega .* timeMesh(:).')' * U );
end