function [m,n] = linTo2D(l,M,N)
m = mod(l-1,N) + 1;
n = floor((l-.1)/(M)) + 1;
end

