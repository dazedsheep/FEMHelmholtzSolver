function u = calcSolution(timeMesh, U, omega)

Nt = length(timeMesh);
Nh = size(U,1);
Ns = size(U,2);

% Harmonic indices
k = (0:Nh-1).';

% Time exponentials (Nt x Nh)
E = exp(1i * omega * timeMesh(:) * k.');

% Zero mode
u = real(E(:,1) * U(1,:));

% Higher harmonics
if Nh > 1
    u = u + 2.*real(E(:,2:end) * U(2:end,:));
end

end
