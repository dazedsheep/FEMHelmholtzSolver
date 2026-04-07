function [x] = projectParameters(x)
etazero = zeros(size(x.refState_1.eta0));
x.refState_1.eta0 = max(etazero,x.refState_1.eta0);
x.refState_2.eta0 = max(etazero,x.refState_2.eta0);
x.refState_3.eta0 = max(etazero,x.refState_3.eta0);
end

