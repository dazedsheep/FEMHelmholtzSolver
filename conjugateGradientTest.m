clear all

% define operators
add = @(a,b) a+b;
minus = @(a,b) a-b;
innerProduct = @(a,b) a'*b;
scalarMul = @(c,a) c.*a;
CGTol = 1e-20;

% our operator A
Amat = [2.04,1;1,2];
b = [1.8;2.3];
A= @(xv) Amat*xv;
M = @(x) x;

% solve min || Ax - b ||^2
[z, iters, res] = conjugateGradient(0, A, b, [0;0], CGTol, 100, innerProduct, add, minus, scalarMul, M);

if sqrt(innerProduct(A(z) - b, A(z) - b)) > 1e-30
    error('Conjugate gradient did not converge.');
else
    fprintf('Conjugate gradient... ok\n');
end