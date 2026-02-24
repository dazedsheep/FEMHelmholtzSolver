function [adjx] = adjointLinearisedForwardOperatorAllAtOnce(elements, timeMeshh, A, x0, x, Gx, Gy)
u = x.modDomain;
% adjoint state in coeffs
adjx.ds = -u .* laplacian(x0.u0, A); 
adjx.db = x0.u0 .* dtt(u, timeMeshh);
adjx.deta = -(x0.u0).^2 .* dtt(u, timeMeshh);

% adjoint w.r.t. du
uh = x0.b0 .* u - 2.* x0.eta0 .* x0.u0 .* u;
dutt = dtt(uh, timeMeshh);

s0laplace = x0.s0 .* laplacian(u, A);
dut = dt(u, timeMeshh);
dutlaplace = laplacian(dut, A);

adjx.du = dutt - s0laplace - dutlaplace;

adjx.dobs = x.obs;

end

