function [P] = calcPressureProfile(omega, T, H, U, cN)
% calc the solution(s) for the multi level harmonic approximation

tind = 0:0.01:0.02;
n_t = size(tind,2);

% compute only for n iterations
iter = 1;

% construct p(t,x)
z = repmat(exp(-1i.*omega.*T.*tind)', 1, size(U,2));
P = zeros(iter, size(z,1), size(z,2));
o = 0;
for j=(min(size(H,1),cN) - iter):min(size(H,1),cN)
    o = o + 1;
    rpC = zeros(size(z,1), size(z,2));
    for k=1:j
        rpC = rpC + squeeze(repmat(H(j,k,:),size(z,1),1)).*repmat(exp(1i.*k.*omega.*tind.*T).', 1, size(U,2));
    end
   
    P(o,:,:) = real(rpC);
end

end

