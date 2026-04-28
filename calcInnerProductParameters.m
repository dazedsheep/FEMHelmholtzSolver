function [val] = calcInnerProductParameters(a,b, elements)

% compute the inner product on the product space and the parameter spaces
[~, a1_eta] = integrate_fun_trimesh(elements.opoints, elements.otri, (a.refState_1.eta0 .* b.refState_1.eta0).');
[~, a1_s]   = integrate_fun_trimesh(elements.opoints, elements.otri, (a.refState_1.s0 .* b.refState_1.s0).');
[~, a1_b]   = integrate_fun_trimesh(elements.opoints, elements.otri, (a.refState_1.b0 .* b.refState_1.b0).');

[~, a2_eta] = integrate_fun_trimesh(elements.opoints, elements.otri, (a.refState_2.eta0 .* b.refState_2.eta0).');
[~, a2_s]   = integrate_fun_trimesh(elements.opoints, elements.otri, (a.refState_2.s0 .* b.refState_2.s0).');
[~, a2_b]   = integrate_fun_trimesh(elements.opoints, elements.otri, (a.refState_2.b0 .* b.refState_2.b0).');

[~, a3_eta] = integrate_fun_trimesh(elements.opoints, elements.otri, (a.refState_3.eta0 .* b.refState_3.eta0).');
[~, a3_s]   = integrate_fun_trimesh(elements.opoints, elements.otri, (a.refState_3.s0 .* b.refState_3.s0).');
[~, a3_b]   = integrate_fun_trimesh(elements.opoints, elements.otri, (a.refState_3.b0 .* b.refState_3.b0).');

val = a1_eta + a1_s + a1_b + a2_eta + a2_s + a2_b + a3_eta + a3_s + a3_b;
end

