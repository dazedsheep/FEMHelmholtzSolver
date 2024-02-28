function [v] = regularizedDirac(gamma,x)
v = zeros(1,size(x,2));
idx = abs(x) <= 2*gamma;

v(idx) = 1/(4*gamma) .* (cos((pi.*x(idx))./(2.*gamma)) + 1);

end

