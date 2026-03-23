function  [kappa] = constructKappaReparameterized(elements, s, b, omega, N)
%% Complex wave number
% also include the zero index
% with space dependent speed of sound and diffusvity

kappa = zeros(size(elements.points,1), N, size(omega,2));
for i = 1:size(omega,2)
    for j = 0:(N-1)
        kappa(:,j+1,i) = omega(i).^2.*b./(s + 1i.*(j).*omega(i));
    end
end

end