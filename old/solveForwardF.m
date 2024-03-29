function [boundaryIndices, elements, p, F] = solveForwardF(elements, c, omega, waveNumber, f, nHarmonics)
%% this function generates boundary data for given point sources in a ball 


% first we calculate p1 using a given excitation

kappa = waveNumber;
beta = 1/c;

n = size(elements.points,1);

h = zeros(n,1);

% not used
g = zeros(1,n);

% consider that 1/c^2 is already included
excitation = 2000*ones(n,1); %constant excitation?

p1 = solveHelmholtzVectorizedTmp(elements, omega, kappa, beta, -excitation, h, g, n);

% figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2), real(p1), 'facecolor', 'interp'); shading interp;
% title("Real part of p_1(x).")
% xlabel('x');
% ylabel('y');

% figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2), imag(p1), 'facecolor', 'interp'); shading interp;
% title("Imag part of the FEM solution using linear lagrange elements.")
% xlabel('x');
% ylabel('y');


% now compute the solutions for the harmonic frequencies

p = zeros(nHarmonics+1, n);
p(1,:) = p1;
F(1,:) = h;
for j = 1:(nHarmonics)
    m = j + 1;
    % not nice
    p_i = zeros(1,n);
    % TODO: vectorize
     for i = 1:j
         p_i = p_i + p(i,:).*p(m-i,:);
     end
    F(m,:) = -f.*m^2.*kappa^2.*p_i';
    p(m,:) = solveHelmholtzVectorizedTmp(elements, m*omega, m*kappa, beta, -F(m,:), h, g, n);
end

%%
boundaryIndices = elements.bedges(:,1);

end

