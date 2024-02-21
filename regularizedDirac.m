function [v] = regularizedDirac(gamma,x)
v = 0;
if abs(x)<=2*gamma
    v = 1/(4*gamma) .* (cos((pi.*x)./(2.*gamma) + 1));
end
end

