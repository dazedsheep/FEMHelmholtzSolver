function u = calcSolution(timeMesh, U, omega)

Nt = length(timeMesh);
Nh = size(U,1);
Ns = size(U,2);

% Harmonic indices
k = (0:Nh-1).';

% Time exponentials (Nt x Nh)
E = exp(1i * omega * timeMesh(:) * k.');

u = real(E(:,1:end) * U(1:end,:));

end
