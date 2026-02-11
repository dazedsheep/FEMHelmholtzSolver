function  [kappa] = constructKappaReparameterized(elements, s, b, omega, N)
%% Complex wave number
% also include the zero index
% with space dependent speed of sound and diffusvity

kappa = zeros(size(elements.points,1), N, size(omega,2));
for i = 1:size(omega,2)
    for j = 2:(N)
        kappa(:,j) = omega(i)./sqrt(s + 1i.*j.*omega(i).*b);
    end
end

end