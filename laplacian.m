function [lapu] = laplacian(u, A)
lapu = zeros(size(u));
for i = 1:size(u,1)
    lapu(i,:)   = (A * u(i,:).');
end
end

