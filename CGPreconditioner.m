function [x] = CGPreconditioner(x,M)
eps = 1e-8;

x.refState_1.s0 = x.refState_1.s0 ./ (M.refState_1.Mlapu0.' + eps);
x.refState_1.b0 = x.refState_1.b0 ./ (M.refState_1.Mu0tt.' + eps);
x.refState_1.eta0 = x.refState_1.eta0 ./ (M.refState_1.Mu0sqtt.' + eps);

x.refState_2.s0 = x.refState_2.s0 ./ (M.refState_2.Mlapu0.' + eps);
x.refState_2.b0 = x.refState_2.b0 ./ (M.refState_2.Mu0tt.' + eps);
x.refState_2.eta0 = x.refState_2.eta0 ./ (M.refState_2.Mu0sqtt.' + eps);

x.refState_3.s0 = x.refState_3.s0 ./ (M.refState_3.Mlapu0.' + eps);
x.refState_3.b0 = x.refState_3.b0 ./ (M.refState_3.Mu0tt.' + eps);
x.refState_3.eta0 = x.refState_3.eta0 ./ (M.refState_3.Mu0sqtt.' + eps);
end

