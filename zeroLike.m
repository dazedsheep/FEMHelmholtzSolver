function z = zeroLike(x)
z.s   = zeros(size(x.s));
z.b   = zeros(size(x.b));
z.eta = zeros(size(x.eta));
z.u   = zeros(size(x.u));
end