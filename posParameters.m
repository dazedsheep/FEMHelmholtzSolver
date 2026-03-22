function [x] = posParameters(x)

x.refState_1.eta0   = max(x.refState_1.eta0,0);
x.refState_1.s0     = max(x.refState_1.s0,0);
x.refState_1.b0     = max(x.refState_1.b0,0);

x.refState_2.eta0   = max(x.refState_2.eta0,0);
x.refState_2.s0     = max(x.refState_2.s0,0);
x.refState_2.b0     = max(x.refState_2.b0,0);

x.refState_3.eta0   = max(x.refState_3.eta0,0);
x.refState_3.s0     = max(x.refState_3.s0,0);
x.refState_3.b0     = max(x.refState_3.b0,0);
end

