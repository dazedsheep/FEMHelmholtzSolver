function[val] = minusParameters(a, b)
% a - b
val = a;
val.refState_1.eta0   = a.refState_1.eta0  - b.refState_1.eta0;
val.refState_1.s0     = a.refState_1.s0    - b.refState_1.s0;
val.refState_1.b0     = a.refState_1.b0    - b.refState_1.b0;

val.refState_2.eta0   = a.refState_2.eta0  - b.refState_2.eta0;
val.refState_2.s0     = a.refState_2.s0    - b.refState_2.s0;
val.refState_2.b0     = a.refState_2.b0    - b.refState_2.b0;

val.refState_3.eta0   = a.refState_3.eta0  - b.refState_3.eta0;
val.refState_3.s0     = a.refState_3.s0    - b.refState_3.s0;
val.refState_3.b0     = a.refState_3.b0    - b.refState_3.b0;

end