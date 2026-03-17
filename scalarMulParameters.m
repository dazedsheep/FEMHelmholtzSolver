function [val] = scalarMulParameters(c,a)
%c ... scalar
%
val = a;
val.refState_1.eta0   = c.*a.refState_1.eta0;  
val.refState_1.s0     = c.*a.refState_1.s0;   
val.refState_1.b0     = c.*a.refState_1.b0; 

val.refState_2.eta0   = c.*a.refState_2.eta0;
val.refState_2.s0     = c.*a.refState_2.s0;   
val.refState_2.b0     = c.*a.refState_2.b0;  

val.refState_3.eta0   = c.*a.refState_3.eta0;
val.refState_3.s0     = c.*a.refState_3.s0;   
val.refState_3.b0     = c.*a.refState_3.b0;    
end

