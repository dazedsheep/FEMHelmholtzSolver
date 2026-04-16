function [b_adj, s_adj, eta_adj] = calcAdjointStates(Uadj, omega, timeMesh, u0Laplace, u0tt, u0Sqtt)

adjointState = calcSolution(timeMesh, Uadj, omega);

s_adj = adjointState.*u0Laplace;
s_adj = trapz(timeMesh.', s_adj,1).';

b_adj = adjointState.*u0tt;
b_adj = -trapz(timeMesh.', b_adj,1).';

eta_adj = adjointState.*u0Sqtt;
eta_adj = trapz(timeMesh.', eta_adj,1).';
end

