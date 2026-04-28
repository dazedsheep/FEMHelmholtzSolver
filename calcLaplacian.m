function DeltaU_time = calcLaplacian(timeMesh, u, kappa, F, omega)

[Np1, n] = size(u);
N = Np1 - 1;
Nt = length(timeMesh);

DeltaU_time = zeros(Nt, n);

for m = 1:N
    
    Delta_um = ...
        - m^2 .* kappa(:,m+1).' .* u(m+1,:) ...
        - F(m+1,:);
    
    phase = exp(1i*m*omega*timeMesh(:));   % Nt x 1
    
    DeltaU_time = DeltaU_time + real( phase * Delta_um );
end

end